#!/usr/bin/env python3

# This file is part of the Mamico project
# Author: Piet Jarmatz
# July 2020
# Helmut Schmidt University, Hamburg. Chair for High Performance Computing
# BSD license, see the copyright notice in Mamico's main folder

import sys
sys.path.append('../../../../build')
# sys.path.append('../../coupling/filtering/filters')
import adios2
import matplotlib.pyplot as mplt
from configparser import ConfigParser
import pandas as pd
import numpy as np
import mamico.tarch.utils
import mamico.coupling
from mamico.coupling.solvers import CouetteSolverInterface
from mamico.coupling.services import MultiMDCellService
import mamico.tarch.configuration
import json
import math
import logging
import coloredlogs

log = logging.getLogger('KVSTest')
logging.getLogger('matplotlib.font_manager').disabled = True


BENCH_BEFORE_RUN = False
RANK = mamico.tarch.utils.initMPI()

# Versatile configurable MPI parallel Kármán vortex street flow test for noise-filtered multi-instance Nie coupling.
# Features:
# -> one-way or two-way coupling
# -> using lbmpy for Lattice Bolzmann flow simulation with CUDA on GPU
# -> using (multi-instance) SimpleMD or synthetic MD data (with Gaussian noise) (TODO)
# -> using filtering subsystem with configurable coupling data analysis or noise filter sequences (TODO)


class KVSTest():
    def __init__(self, cfg):
        self.cfg = cfg
        self.rank = RANK
        if self.rank == 0:
            log.info("Created KVSTest ...")

    def run(self):
        self.parseXMLConfigurations()
        self.initSolvers()
        for cycle in range(self.cfg.getint("coupling", "couplingCycles")):
            self.runOneCouplingCycle(cycle)
        pd.DataFrame(self.velLB).to_csv("lbm.csv", sep=";", header=False)

    def runOneCouplingCycle(self, cycle):
        self.advanceMacro(cycle)
        self.advanceMicro(cycle)
        self.twoWayCoupling(cycle)
        if self.rank == 0 and (cycle % 20 == 0 or cycle < 10):
            log.info("Finish coupling cycle " + str(cycle))

    def parseXMLConfigurations(self):
        self.simpleMDConfig = mamico.tarch.configuration.parseMolecularDynamicsConfiguration(
            "kvs.xml", "molecular-dynamics")
        if not self.simpleMDConfig.isValid():
            log.error("Invalid SimpleMD config!")
            sys.exit(1)
        self.mamicoConfig = mamico.tarch.configuration.parseMaMiCoConfiguration(
            "kvs.xml", "mamico")
        if not self.mamicoConfig.isValid():
            log.error("Invalid MaMiCo config!")
            sys.exit(1)

    def initSolvers(self):
        if self.rank == 0:
            self.macroscopicSolver = LBSolver(self.cfg)

            # Compute time step size / coupling scaling
            # viscosity must exactly match on MD and LB side

            # to track time passed so far, in Mamico Units
            self.t = 0
            # kinematic viscosity of fluid, in Mamico Units
            kinVisc = self.cfg.getfloat(
                "macroscopic-solver", "viscosity") / self.cfg.getfloat("microscopic-solver", "density")
            # dx LB in Mamico Units
            self.dx = self.mamicoConfig.getMacroscopicCellConfiguration(
            ).getMacroscopicCellSize()[0]
            # time intervall of one coupling cycle, in Mamico Units
            self.dt = self.simpleMDConfig.getSimulationConfiguration().getDt() * \
                self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
            # duration of one LB timestep, in Mamico Units
            self.dt_LB = self.dx*self.dx * \
                (1/self.macroscopicSolver.omega-0.5)/(3*kinVisc)

            log.info("MD timesteps per coupling cycle = "
                     + str(self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()))
            log.info("LB timesteps per coupling cycle = "
                     + str(self.dt / self.dt_LB))

            isteps = self.cfg.getint("macroscopic-solver", "init-timesteps")
            log.info("Running " + str(isteps)
                     + " LB initialisation timesteps ...")

            if BENCH_BEFORE_RUN == True:
                timeguess = isteps * self.macroscopicSolver.domain_size[0] * \
                    self.macroscopicSolver.domain_size[1] * \
                    self.macroscopicSolver.domain_size[2] / \
                    (self.macroscopicSolver.mlups * 1e6)
                log.info("(Estimated runtime = " + str(timeguess) + " seconds)")

            self.macroscopicSolver.advance(isteps)
            log.info("Finished LB init-timesteps.")

        numMD = self.cfg.getint("microscopic-solver", "number-md-simulations")

        self.multiMDService = mamico.tarch.utils.MultiMDService(numberProcesses=self.simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),
                                                                totalNumberMDSimulations=numMD)

        self.localMDInstances = self.multiMDService.getLocalNumberOfMDSimulations()

        if self.rank == 0:
            log.info("totalNumberMDSimulations = " + str(numMD))
            log.info("localMDInstances on rank 0 = "
                     + str(self.localMDInstances))

        self.simpleMD = [mamico.coupling.getMDSimulation(self.simpleMDConfig, self.mamicoConfig,
                                                         self.multiMDService.getLocalCommunicator()) for i in range(self.localMDInstances)]

        for i in range(self.localMDInstances):
            self.simpleMD[i].init(
                self.multiMDService, self.multiMDService.getGlobalNumberOfLocalMDSimulation(i))

        equSteps = self.cfg.getint("microscopic-solver", "equilibration-steps")
        for i in range(self.localMDInstances):
            self.simpleMD[i].switchOffCoupling()
            self.simpleMD[i].simulateTimesteps(equSteps, 0)
            self.simpleMD[i].switchOnCoupling()
        self.mdStepCounter = equSteps

        #Warning: If instances < ranks, this is empty for rank 0. In that case, the MultiMDCellService intialisation below segfaults. FIXME
        self.mdSolverInterface = [mamico.coupling.getMDSolverInterface(self.simpleMDConfig, self.mamicoConfig,
                                                                       self.simpleMD[i]) for i in range(self.localMDInstances)]

        self.macroscopicSolverInterface = CouetteSolverInterface(self.getGlobalNumberMacroscopicCells(),
                                                                 self.mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap())
        mamico.coupling.setMacroscopicSolverInterface(
            self.macroscopicSolverInterface)

        self.multiMDCellService = MultiMDCellService(self.mdSolverInterface, self.macroscopicSolverInterface,
                                                     self.simpleMDConfig, self.rank, self.cfg.getint(
                                                         "microscopic-solver", "number-md-simulations"),
                                                     self.mamicoConfig, "kvs.xml", self.multiMDService)

        mamico.tarch.utils.initIndexing(
            self.simpleMDConfig, self.mamicoConfig, self.macroscopicSolverInterface, self.rank)

        self.multiMDCellService.constructFilterPipelines()

        for i in range(self.localMDInstances):
            self.simpleMD[i].setMacroscopicCellService(
                self.multiMDCellService.getMacroscopicCellService(i))
            self.multiMDCellService.getMacroscopicCellService(i).computeAndStoreTemperature(
                self.cfg.getfloat("microscopic-solver", "temperature"))

        from scipy.ndimage import gaussian_filter, median_filter
        #Add Gauss filter

        def gauss_sca05(data):
            print("Applying gaussian filter to a scalar property. sigma = 0.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(0.5, 0.5, 0.5))

        def gauss_vec05(data):
            print("Applying gaussian filter to a 3d property. sigma = 0.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(0.5, 0.5, 0.5, 0))

        def gauss_sca1(data):
            print("Applying gaussian filter to a scalar property. sigma = 1.")
            return gaussian_filter(data, truncate=1.0, sigma=(1, 1, 1), mode="mirror")

        def gauss_vec1(data):
            print("Applying gaussian filter to a 3d property. sigma = 1.")
            return gaussian_filter(data, truncate=1.0, sigma=(1, 1, 1, 0), mode="mirror")

        def gauss_sca15(data):
            print("Applying gaussian filter to a scalar property. sigma = 1.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(1.5, 1.5, 1.5))

        def gauss_vec15(data):
            print("Applying gaussian filter to a 3d property. sigma = 1.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(1.5, 1.5, 1.5, 0))

        def gauss_sca2(data):
            print("Applying gaussian filter to a scalar property. sigma = 2.")
            return gaussian_filter(data, truncate=1.0, sigma=(2, 2, 2))

        def gauss_vec2(data):
            print("Applying gaussian filter to a 3d property. sigma = 2.")
            return gaussian_filter(data, truncate=1.0, sigma=(2, 2, 2, 0))

        def gauss_sca25(data):
            print("Applying gaussian filter to a scalar property. sigma = 2.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(2.5, 2.5, 2.5))

        def gauss_vec25(data):
            print("Applying gaussian filter to a 3d property. sigma = 2.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(2.5, 2.5, 2.5, 0))

        def gauss_sca3(data):
            print("Applying gaussian filter to a scalar property. sigma = 3.")
            return gaussian_filter(data, truncate=1.0, sigma=(3, 3, 3))

        def gauss_vec3(data):
            print("Applying gaussian filter to a 3d property. sigma = 3.")
            return gaussian_filter(data, truncate=1.0, sigma=(3, 3, 3, 0))

        def gauss_sca35(data):
            print("Applying gaussian filter to a scalar property. sigma = 3.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(3.5, 3.5, 3.5))

        def gauss_vec35(data):
            print("Applying gaussian filter to a 3d property. sigma = 3.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(3.5, 3.5, 3.5, 0))

        def gauss_sca4(data):
            print("Applying gaussian filter to a scalar property. sigma = 4.")
            return gaussian_filter(data, truncate=1.0, sigma=(4, 4, 4))

        def gauss_vec4(data):
            print("Applying gaussian filter to a 3d property. sigma = 4.")
            return gaussian_filter(data, truncate=1.0, sigma=(4, 4, 4, 0))

        def gauss_sca45(data):
            print("Applying gaussian filter to a scalar property. sigma = 4.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(4.5, 4.5, 4.5))

        def gauss_vec45(data):
            print("Applying gaussian filter to a 3d property. sigma = 4.5.")
            return gaussian_filter(data, truncate=1.0, sigma=(4.5, 4.5, 4.5, 0))

        def gauss_sca5(data):
            print("Applying gaussian filter to a scalar property. sigma = 5.")
            return gaussian_filter(data, truncate=1.0, sigma=(5, 5, 5))

        def gauss_vec5(data):
            print("Applying gaussian filter to a 3d property. sigma = 5.")
            return gaussian_filter(data, truncate=1.0, sigma=(5, 5, 5, 0))

        def gauss_x_sca1(data):
            print(
                "Applying gaussian filter to a scalar property. Filtering on X-axis only. sigma = 1.")
            return gaussian_filter(data, truncate=1.0, sigma=(1, 0, 0), mode="mirror")

        def gauss_x_vec1(data):
            print(
                "Applying gaussian filter to a 3d property. Filtering on X-axis only. sigma = 1.")
            return gaussian_filter(data, truncate=1.0, sigma=(1, 0, 0, 0), mode="mirror")

        mcs = self.multiMDCellService.getMacroscopicCellService(0)

        #fff testing
        #mcs.addFilterToSequence(filter_sequence="gauss-wtf", filter_index=0, scalar_filter_func = gauss_sca1, vector_filter_func=gauss_vec1)
        #mcs.addFilterToSequence(filter_sequence="gauss", filter_index=0, scalar_filter_func = gauss_sca1, vector_filter_func=gauss_vec1)

        #gauss testing
        #mcs.addFilterToSequence(filter_sequence="gaussX-python", filter_index=0, scalar_filter_func = gauss_x_sca1, vector_filter_func=gauss_x_vec1)
        #mcs.addFilterToSequence(filter_sequence="gaussXYZ-python", filter_index=0, scalar_filter_func = gauss_sca1, vector_filter_func=gauss_vec1)

        #param testing
        #mcs.addFilterToSequence(filter_sequence="gauss-05", filter_index=0, scalar_filter_func = gauss_sca05, vector_filter_func=gauss_vec05)
        #mcs.addFilterToSequence(filter_sequence="gauss-1", filter_index=0, scalar_filter_func = gauss_sca1, vector_filter_func=gauss_vec1)
        #mcs.addFilterToSequence(filter_sequence="gauss-15", filter_index=0, scalar_filter_func = gauss_sca15, vector_filter_func=gauss_vec15)
        #mcs.addFilterToSequence(filter_sequence="gauss-2", filter_index=0, scalar_filter_func = gauss_sca2, vector_filter_func=gauss_vec2)
        #mcs.addFilterToSequence(filter_sequence="gauss-25", filter_index=0, scalar_filter_func = gauss_sca25, vector_filter_func=gauss_vec25)
        #mcs.addFilterToSequence(filter_sequence="gauss-3", filter_index=0, scalar_filter_func = gauss_sca3, vector_filter_func=gauss_vec3)
        #mcs.addFilterToSequence(filter_sequence="gauss-35", filter_index=0, scalar_filter_func = gauss_sca35, vector_filter_func=gauss_vec35)
        #mcs.addFilterToSequence(filter_sequence="gauss-4", filter_index=0, scalar_filter_func = gauss_sca4, vector_filter_func=gauss_vec4)
        #mcs.addFilterToSequence(filter_sequence="gauss-45", filter_index=0, scalar_filter_func = gauss_sca45, vector_filter_func=gauss_vec45)
        #mcs.addFilterToSequence(filter_sequence="gauss-5", filter_index=0, scalar_filter_func = gauss_sca5, vector_filter_func=gauss_vec5)

        self.buf = mamico.coupling.Buffer(self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),
                                          self.macroscopicSolverInterface, self.rank, self.mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap())

        self.csv = self.cfg.getint("coupling", "csv-every-timestep")
        self.png = self.cfg.getint("coupling", "png-every-timestep")
        self.adios2 = self.cfg.getint("coupling", "adios2-every-timestep")

        # buffer for evaluation plot
        if self.rank == 0:
            self.velLB = np.zeros(
                (self.cfg.getint("coupling", "couplingCycles"), 2))
            if self.adios2 > 0:
                self.adiosfile = adios2.open("kvstest_volume.bp", "w")
                timefactor = self.adios2 * self.dt/(self.simpleMDConfig.getADIOS2Configuration(
                ).getWriteEveryTimestep() * self.simpleMDConfig.getSimulationConfiguration().getDt())
                print("timefactor:", timefactor)
                self.adiosfile.write_attribute('timefactor', str(timefactor))
        if self.rank == 0:
            log.info("Finished initSolvers")  # after ? ms

    def shutdown(self):
        if self.rank == 0:
            log.info("Finished " + str(self.mdStepCounter) + " MD timesteps")
            if self.adios2 > 0:
                self.adiosfile.close()
        #Analyse data gathered by StrouhalPython filter
        #self.sf.calculateStrouhalNumber()

        for i in range(self.localMDInstances):
            self.simpleMD[i].shutdown()
        # TODO something else to do here??

    def advanceMacro(self, cycle):
        if self.rank == 0:
            to_advance = (cycle+1)*self.dt - self.t
            steps = int(round(to_advance / self.dt_LB))
            log.debug("Advancing " + str(steps) + " LB timesteps")
            self.macroscopicSolver.advance(steps)
            self.t = self.t + steps * self.dt_LB

            cellmass = (self.cfg.getfloat("microscopic-solver", "density")
                        * np.prod(self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()))
            numcells = self.getGlobalNumberMacroscopicCells()
            mdpos = json.loads(self.cfg.get("domain", "md-pos"))
            # Convert center of MD domain in SI units to offset of MD domain as cell index
            mdpos = [int(mdpos[d]*self.macroscopicSolver.cpm - numcells[d]/2)
                     for d in range(3)]

            self.buf.store2send(cellmass,
                                self.macroscopicSolver.scen.velocity[
                                    mdpos[0]:mdpos[0]+numcells[0],
                                    mdpos[1]:mdpos[1]+numcells[1],
                                    mdpos[2]:mdpos[2]+numcells[2],
                                    :].data * (self.dx / self.dt_LB),
                                self.macroscopicSolver.scen.density[
                                    mdpos[0]:mdpos[0]+numcells[0],
                                    mdpos[1]:mdpos[1]+numcells[1],
                                    mdpos[2]:mdpos[2]+numcells[2]].data
                                )

            if self.png > 0 and (cycle+1) % self.png == 0:
                filename = "kvstest_" + str(cycle+1) + ".png"
                self.macroscopicSolver.plot(filename)

            if self.adios2 > 0 and (cycle+1) % self.adios2 == 0:
               to_write = np.ascontiguousarray(
                   self.macroscopicSolver.scen.velocity[:, :, :, :].data, dtype=np.float32)
               log.info("writing to adios2 " + str(type(to_write)) + " with shape "
                        + str(to_write.shape) + " and dtype " + str(to_write.dtype))
               shape = np.array(self.macroscopicSolver.scen.velocity[:, :, :, 0].data.shape) * \
                   self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()
               offset = [0., 0., 0.]
               print(shape)
               print(mdpos)
               print(self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize())
               offset[0] = -((mdpos[0]/to_write.shape[0] * shape[0])-(
                   0.5 * self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()[0]))
               offset[1] = -((mdpos[1]/to_write.shape[1] * shape[1])-(
                   0.5 * self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()[1]))
               offset[2] = -((mdpos[2]/to_write.shape[2] * shape[2])-(
                   0.5 * self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()[2]))
               print(offset)
               gb_to_write = np.array(
                   [offset[0], offset[1], offset[2], offset[0] + shape[0], offset[1] + shape[1], offset[2] + shape[2]])
               print(gb_to_write)
               self.adiosfile.write("global_box", gb_to_write, gb_to_write.shape, np.zeros_like(
                   gb_to_write.shape), gb_to_write.shape)
               self.adiosfile.write("velocity", to_write, to_write.shape, np.zeros_like(
                   to_write.shape), to_write.shape)
               self.adiosfile.end_step()

        if self.cfg.getboolean("coupling", "send-from-macro-to-md"):
            self.multiMDCellService.sendFromMacro2MD(self.buf)

    def advanceMicro(self, cycle):
        numT = self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
        for i in range(self.localMDInstances):
            mamico.coupling.setMacroscopicCellService(
                self.multiMDCellService.getMacroscopicCellService(i))
            mamico.coupling.setMDSolverInterface(self.mdSolverInterface[i])
            self.simpleMD[i].simulateTimesteps(numT, self.mdStepCounter)
            self.multiMDCellService.getMacroscopicCellService(
                i).plotEveryMacroscopicTimestep(cycle)
        self.mdStepCounter = self.mdStepCounter + numT

        if self.cfg.getboolean("coupling", "send-from-md-to-macro"):
            self.multiMDCellService.sendFromMD2Macro(self.buf)

    def twoWayCoupling(self, cycle):
        if self.cfg.getboolean("coupling", "two-way-coupling") and self.rank == 0:
            self.macroscopicSolver.setMDBoundaryValues(self.buf)  # TODO

        if self.csv > 0 and (cycle+1) % self.csv == 0 and self.rank == 0:
            filename = "KVSMD2Macro_" + str(cycle+1) + ".csv"
            self.buf.recv2CSV(filename)

        if self.rank == 0:
            # extract data for eval plots
            numcells = self.getGlobalNumberMacroscopicCells()
            mdpos = json.loads(self.cfg.get("domain", "md-pos"))
            # Convert center of MD domain in SI units to offset of MD domain as cell index
            mdpos = [int(mdpos[d]*self.macroscopicSolver.cpm - numcells[d]/2)
                     for d in range(3)]
            for dir in range(2):
                self.velLB[cycle, dir] = self.macroscopicSolver.scen.velocity[
					#TODO: exact cell index in md2macro
                    mdpos[0]+6,
                    mdpos[1]+6,
                    mdpos[2]+6,
                    dir].data * (self.dx / self.dt_LB)

    def __del__(self):
        self.shutdown()
        mamico.tarch.utils.finalizeMPI()
        if self.rank == 0:
            log.info("Finished KVSTest shutdown")

    def getGlobalNumberMacroscopicCells(self):
        domainSize = self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()
        dx = self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()
        return [math.floor(domainSize[d]/dx[d]+0.5) for d in range(3)]


if RANK == 0:    # This fixes last_config.json-Error
    from lbmpy.session import *
    from lbmpy.parameterization import Scaling


lb_log = logging.getLogger('LBSolver')


class LBSolver():
    def __init__(self, cfg):
        self.cfg = cfg
        self.cpm = cpm = cfg.getint("macroscopic-solver", "cells-per-meter")
        self.domain_size = (int(2.5*cpm), int(0.41*cpm), int(0.41*cpm))
        self.vis = 1e-3
        self.scaling = Scaling(physical_length=0.1, physical_velocity=2.25, kinematic_viscosity=self.vis,
                               cells_per_length=0.1*cpm)
        self.omega = cfg.getfloat("macroscopic-solver", "omega")
        if self.omega > 1.92:
            lb_log.warning("High omega, LB may be unstable!")
        self.scaling_result = self.scaling.diffusive_scaling(self.omega)
        lb_log.info("Scaling LB to SI units:")
        lb_log.info(self.scaling_result)
        lb_log.info("dx=" + str(self.scaling.dx))
        if self.scaling_result.lattice_velocity > 0.05:
            lb_log.warning("High lattice velocity, LB may be unstable!")
        if self.scaling_result.lattice_velocity > 0.57:
            lb_log.error("Supersonic lattice velocity!")

        self.timesteps_finished = 0
        self.setup_scenario()
        lb_log.info("Successfully created scenario")
        lb_log.info("Domain size = " + str(self.domain_size))
        lb_log.info("Total number of cells = "
                    + str(self.domain_size[0]*self.domain_size[1]*self.domain_size[2]))

        if BENCH_BEFORE_RUN == True:
            lb_log.info("Running benchmark ...")
            self.mlups = self.scen.benchmark()
       	    lb_log.info("Benchmark result = " + str(self.mlups) + " MLUPS")

        self.cD_max = 0
        self.cL_max = 0

    def setup_scenario(self):
        try:  # for lbmpy version >= 0.4.0
            from pystencils import Target
            if self.cfg.get("macroscopic-solver", "optimization-target") == "cpu":
                optTarget = Target.CPU
            elif self.cfg.get("macroscopic-solver", "optimization-target") == "gpu":
                optTarget = Target.GPU
        except ImportError:   # for lbmpy version <= 0.3.4
            optTarget = self.cfg.get(
                "macroscopic-solver", "optimization-target")
        self.scen = LatticeBoltzmannStep(domain_size=self.domain_size, method='trt', stencil='D3Q19',
                                         relaxation_rate=self.omega, periodicity=(
                                             True, False, False),
                                         optimization={'target': optTarget,
                                                       'gpu_indexing': 'line',
                                                       'double_precision': self.cfg.get("macroscopic-solver", "double-precision")})

        def obstacle(x, y, z):
            return (x > 0.45*self.cpm) & (x < 0.55*self.cpm) & (y > 0.15*self.cpm) & (y < 0.25*self.cpm)
        wall = NoSlip()
        self.scen.boundary_handling.set_boundary(wall, mask_callback=obstacle)

        def velocity_info_callback(boundary_data, **_):
            boundary_data['vel_1'] = 0
            boundary_data['vel_2'] = 0
            u_max = self.scaling_result.lattice_velocity
            y, z = boundary_data.link_positions(
                1), boundary_data.link_positions(2)
            H = self.domain_size[1]
            boundary_data['vel_0'] = 16 * u_max * \
                y * z * (H-y) * (H-z) / (H*H*H*H)
        inflow = UBB(velocity_info_callback, dim=self.scen.method.dim)
        self.scen.boundary_handling.set_boundary(inflow, make_slice[0, :, :])

        outflow = FixedDensity(1.0)
        self.scen.boundary_handling.set_boundary(outflow, make_slice[-1, :, :])

    def advance(self, timesteps):
        vtk_every_ts = 100000
       	ts_goal = self.timesteps_finished + timesteps
        while self.timesteps_finished < ts_goal:
            if ts_goal - self.timesteps_finished < vtk_every_ts:
                steps = ts_goal - self.timesteps_finished
                self.scen.run(steps)
            else:
                steps = vtk_every_ts
                self.scen.run(steps)
                self.scen.write_vtk()
            self.timesteps_finished = self.timesteps_finished + steps
        self.compute_drag_lift()

    def compute_drag_lift(self):
        #scaling
        dx = self.scaling.dx
        dt = self.scaling_result.dt

        S = 0.1 * 0.41     # surface of one cylinder face in m^2

        # compute cell indices around cylinder obstacle
        # TODO only tested for cpm = 312 cells per meter
        left = int(0.45*self.cpm) - 1  # X: left = fluid,   left+1 = obstacle
        right = int(0.55*self.cpm) + 1  # X: right = fluid,  right-1 = obstacle
        bottom = int(0.15*self.cpm)    # Y: bottom = fluid, bottom+1 = obstacle
        top = int(0.25*self.cpm)       # Y: top = fluid,    top-1 = obstacle

        # compute drag
        kvs = self.scen
        # compute advective drag force
        # pressure = density * c_s^2 = density / 3
        p_left = np.mean(kvs.density[left, bottom+1:top-1, :]) / 3
        p_right = np.mean(kvs.density[right, bottom+1:top-1, :]) / 3
        # convert from lattice to SI units
        p_left = p_left * dx*dx/(dt*dt)
        p_right = p_right * dx*dx/(dt*dt)
        FD = (p_left - p_right) * S
        # compute viscous drag force
        # ignore 3 cell layers, because simple wall boundary condition resets all PDFs to equilibrium => wrong results very close to boundary
        # use second order forward finite difference to get first derivative in normal direction of tangential velocity
        vel_grad_bottom = \
            -3/2 * np.mean(kvs.velocity[left+1:right-1, bottom-3, :, 0]) + \
            2 * np.mean(kvs.velocity[left+1:right-1, bottom-4, :, 0]) + \
            -1/2 * np.mean(kvs.velocity[left+1:right-1, bottom-5, :, 0])
        vel_grad_top = \
            -3/2 * np.mean(kvs.velocity[left+1:right-1, top+3, :, 0]) + \
            2 * np.mean(kvs.velocity[left+1:right-1, top+4, :, 0]) + \
            -1/2 * np.mean(kvs.velocity[left+1:right-1, top+5, :, 0])
        lb_log.debug("vel_grad_bottom = " + str(vel_grad_bottom))
        lb_log.debug("vel_grad_top = " + str(vel_grad_top))
        #convert to SI, integrate over surface, add to force
        FD = FD + self.vis * (vel_grad_bottom + vel_grad_top) / dt * S

        # compute lift
        # compute advective lift force
        p_bottom = np.mean(kvs.density[left+1:right-1, bottom, :]) / 3
        p_top = np.mean(kvs.density[left+1:right-1, top, :]) / 3
        p_bottom = p_bottom * dx*dx/(dt*dt)
        p_top = p_top * dx*dx/(dt*dt)
        FL = (p_bottom - p_top) * S
        # compute viscous lift force
        vel_grad_left = \
            -3/2 * np.mean(kvs.velocity[left-3, bottom+1:top-1, :, 1]) + \
            2 * np.mean(kvs.velocity[left-4, bottom+1:top-1, :, 1]) + \
            -1/2 * np.mean(kvs.velocity[left-5, bottom+1:top-1, :, 1])
        vel_grad_right = \
            -3/2 * np.mean(kvs.velocity[right+3, bottom+1:top-1, :, 1]) + \
            2 * np.mean(kvs.velocity[right+4, bottom+1:top-1, :, 1]) + \
            -1/2 * np.mean(kvs.velocity[right+5, bottom+1:top-1, :, 1])
        FL = FL + self.vis * (vel_grad_left + vel_grad_right) / dt * S

        # Convert forces in Newton to (dimensionless) coefficients
        cD = 2 * FD / (1*1*1*S)
        if(abs(cD) > abs(self.cD_max)):
            self.cD_max = cD

        cL = 2 * FL / (1*1*1*S)
        if(abs(cL) > abs(self.cL_max)):
            self.cL_max = cL

    def __del__(self):
        lb_log.info("Finished " + str(self.timesteps_finished)
                    + " LB timesteps")
        lb_log.info("cD_max = " + str(self.cD_max))
        lb_log.info("cL_max = " + str(self.cL_max))

    # write velocity slice to png file
    def plot(self, filename):
        plt.vector_field_magnitude(
            self.scen.velocity[:, :, int(self.domain_size[2]//2), 0:2])
        plt.subplots_adjust(left=0, bottom=0, right=1,
                            top=1, wspace=None, hspace=None)
        plt.axis('off')
        plt.savefig(filename)

# Read config file, create a KVSTest instance and run it


def main():
    cfg = ConfigParser()
    cfg.read("kvstest.ini")

    # console mode, colored log to stdout
    if len(sys.argv) == 1:
        coloredlogs.install(
            fmt='%(asctime)s.%(msecs)03d %(name)s %(levelname)s %(message)s', level='DEBUG')
        log.setLevel(level=logging.INFO)
        lb_log.setLevel(level=logging.INFO)
    # job mode, log to file
    elif len(sys.argv) == 2:
        logging.basicConfig(filename=sys.argv[1], level=logging.INFO)
    else:
        print("Usage:")
        print("[1] " + str(sys.argv[0]))
        print("[2] " + str(sys.argv[0]) + " logfile")
        sys.exit(1)

    test = KVSTest(cfg)
    test.run()


# only if this file is executed directly ie. not imported as a module,
# then actually run the simulation
if __name__ == '__main__':
    main()
