#!/usr/bin/env python3

# This file is part of the Mamico project
# Author: Piet Jarmatz
# July 2020
# Helmut Schmidt University, Hamburg. Chair for High Performance Computing
# BSD license, see the copyright notice in Mamico's main folder

import sys
sys.path.append('../../python-binding')
import coloredlogs, logging
import math
import json
import mamico.tarch.configuration
from mamico.coupling.services import MultiMDCellService
from mamico.coupling.solvers import CouetteSolverInterface
import mamico.coupling
import mamico.tarch.utils
import numpy as np
from configparser import ConfigParser

log = logging.getLogger('KVSTest')

# Versatile configurable MPI parallel Kármán vortex street flow test for noise-filtered multi-instance Nie coupling.
# Features:
# -> one-way or two-way coupling
# -> using lbmpy for Lattice Bolzmann flow simulation with CUDA on GPU
# -> using (multi-instance) SimpleMD or synthetic MD data (with Gaussian noise) (TODO)
# -> using filtering subsystem with configurable coupling data analysis or noise filter sequences (TODO)
class KVSTest():
    def __init__(self, cfg):
        self.cfg = cfg
        self.rank = mamico.tarch.utils.initMPI()
        if self.rank==0:
            log.info("Created KVSTest ...")

    def run(self):
        self.parseXMLConfigurations()
        self.initSolvers()
        for cycle in range(self.cfg.getint("coupling", "couplingCycles")):
            self.runOneCouplingCycle(cycle)

    def runOneCouplingCycle(self, cycle):
        self.advanceMacro(cycle)
        self.advanceMicro(cycle)
        self.twoWayCoupling(cycle)
        if self.rank==0 and (cycle%20==0 or cycle<10):
            log.info("Finish coupling cycle " +str(cycle))
    
    def parseXMLConfigurations(self):
        self.simpleMDConfig = mamico.tarch.configuration.parseMolecularDynamicsConfiguration("couette_simplemd.xml","molecular-dynamics")
        if not self.simpleMDConfig.isValid():
            log.error("Invalid SimpleMD config!")
            sys.exit(1)
        self.mamicoConfig = mamico.tarch.configuration.parseMaMiCoConfiguration("couette_mamico.xml","mamico")
        if not self.mamicoConfig.isValid():
            log.error("Invalid MaMiCo config!")
            sys.exit(1)

    def initSolvers(self):
        if self.rank==0:
            self.macroscopicSolver = LBSolver(self.cfg)

            # Compute time step size / coupling scaling 
            # viscosity must exactly match on MD and LB side

            # to track time passed so far, in Mamico Units
            self.t = 0 
            # kinematic viscosity of fluid, in Mamico Units
            kinVisc = self.cfg.getfloat("macroscopic-solver", "viscosity") / self.cfg.getfloat("microscopic-solver", "density");
            # dx LB in Mamico Units
            self.dx = self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0]
            # time intervall of one coupling cycle, in Mamico Units
            self.dt = self.simpleMDConfig.getSimulationConfiguration().getDt() * \
                self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
            # duration of one LB timestep, in Mamico Units
            self.dt_LB = self.dx*self.dx*(1/self.macroscopicSolver.omega-0.5)/(3*kinVisc)

            log.info("MD timesteps per coupling cycle = " + 
                str(self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()))
            log.info("LB timesteps per coupling cycle = " + str(self.dt / self.dt_LB))

        numMD = self.cfg.getint("microscopic-solver","number-md-simulations")

        self.multiMDService = mamico.tarch.utils.MultiMDService(numberProcesses = self.simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),
            totalNumberMDSimulations = numMD)

        self.localMDInstances = self.multiMDService.getLocalNumberOfMDSimulations()
        if self.rank==0:
            log.info("totalNumberMDSimulations = " + str(numMD))
            log.info("localMDInstances on rank 0 = " + str(self.localMDInstances))

        self.simpleMD = [mamico.coupling.getMDSimulation(self.simpleMDConfig,self.mamicoConfig, 
            self.multiMDService.getLocalCommunicator()) for i in range(self.localMDInstances)]

        for i in range(self.localMDInstances):
            self.simpleMD[i].init(self.multiMDService, self.multiMDService.getGlobalNumberOfLocalMDSimulation(i))
        
        equSteps = self.cfg.getint("microscopic-solver", "equilibration-steps")
        for i in range(self.localMDInstances):
            self.simpleMD[i].switchOffCoupling()
            self.simpleMD[i].simulateTimesteps(equSteps,0)
            self.simpleMD[i].switchOnCoupling()
        self.mdStepCounter = equSteps
    
        self.mdSolverInterface = [mamico.coupling.getMDSolverInterface(self.simpleMDConfig, self.mamicoConfig,
            self.simpleMD[i]) for i in range(self.localMDInstances)]

        self.macroscopicSolverInterface = CouetteSolverInterface(self.getGlobalNumberMacroscopicCells(),
            self.mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap())
        mamico.coupling.setMacroscopicSolverInterface(self.macroscopicSolverInterface)
    
        self.multiMDCellService = MultiMDCellService(self.mdSolverInterface, self.macroscopicSolverInterface, 
            self.simpleMDConfig, self.rank, self.cfg.getint("microscopic-solver","number-md-simulations"),
            self.mamicoConfig,self.multiMDService)

        for i in range(self.localMDInstances):
            self.simpleMD[i].setMacroscopicCellService(self.multiMDCellService.getMacroscopicCellService(i))
            self.multiMDCellService.getMacroscopicCellService(i).computeAndStoreTemperature(
                self.cfg.getfloat("microscopic-solver", "temperature"))
      
        self.buf = mamico.coupling.Buffer(self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),
            self.macroscopicSolverInterface, self.rank)

        self.csv = self.cfg.getint("coupling", "csv-every-timestep")
        self.png = self.cfg.getint("coupling", "png-every-timestep")
    
        if self.rank==0:
            log.info("Finished initSolvers") # after ? ms

    def shutdown(self):
        if self.rank==0:
            log.info("Finished " + str(self.mdStepCounter) + " MD timesteps")
        for i in range(self.localMDInstances):
            self.simpleMD[i].shutdown()
        # TODO something else to do here??

    def advanceMacro(self,cycle):
        if self.rank==0:
            to_advance = (cycle+1)*self.dt - self.t
            steps = int(round(to_advance / self.dt_LB))
            log.debug("Advancing " + str(steps) + " LB timesteps")
            self.macroscopicSolver.advance(steps) 
            self.t = self.t + steps * self.dt_LB

            cellmass = ( self.cfg.getfloat("microscopic-solver", "density")
                * np.prod(self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()))
            numcells = self.getGlobalNumberMacroscopicCells()
            mdpos = json.loads(self.cfg.get("domain", "md-pos"))
            # Convert center of MD domain in SI units to offset of MD domain as cell index
            mdpos = [int(mdpos[d]*self.macroscopicSolver.cpm - numcells[d]/2) 
                for d in range(3)]

            self.buf.store2send(cellmass, self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),
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

            if self.png > 0 and (cycle+1)%self.png == 0:
                filename = "kvstest_" + str(cycle+1) + ".png"
                self.macroscopicSolver.plot(filename)

        if self.cfg.getboolean("coupling", "send-from-macro-to-md"):
            self.multiMDCellService.sendFromMacro2MD(self.buf)
    
    def advanceMicro(self, cycle):
        numT = self.simpleMDConfig.getSimulationConfiguration().getNumberOfTimesteps()
        for i in range(self.localMDInstances):
            mamico.coupling.setMacroscopicCellService(self.multiMDCellService.getMacroscopicCellService(i))
            mamico.coupling.setMDSolverInterface(self.mdSolverInterface[i])
            self.simpleMD[i].simulateTimesteps(numT,self.mdStepCounter)
            self.multiMDCellService.getMacroscopicCellService(i).plotEveryMacroscopicTimestep(cycle)
        self.mdStepCounter = self.mdStepCounter + numT

        if self.cfg.getboolean("coupling", "send-from-md-to-macro"):
            self.multiMDCellService.sendFromMD2Macro(self.buf)

    def twoWayCoupling(self, cycle):
        idcv = self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion()
        if self.cfg.getboolean("coupling", "two-way-coupling") and self.rank==0:
            self.macroscopicSolver.setMDBoundaryValues(self.buf, idcv) # TODO
        
        if self.csv > 0 and (cycle+1)%self.csv == 0 and self.rank==0:
            filename = "KVSMD2Macro_" + str(cycle+1) + ".csv"
            self.buf.recv2CSV(idcv,filename)

    def __del__(self):
        self.shutdown()
        mamico.tarch.utils.finalizeMPI()
        if self.rank==0:
            log.info("Finished KVSTest shutdown")

    def getGlobalNumberMacroscopicCells(self):
        domainSize = self.simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()
        dx = self.mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()
        return [math.floor(domainSize[d]/dx[d]+0.5) for d in range(3)]

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
        lb_log.info("Total number of cells = " + str(self.domain_size[0]*self.domain_size[1]*self.domain_size[2]))
        lb_log.info("Running benchmark ...")
        lb_log.info("Benchmark result = " + str(self.scen.benchmark()) + " MLUPS")

        self.cD_max = 0
        self.cL_max = 0

    def setup_scenario(self):
        self.scen = LatticeBoltzmannStep(domain_size=self.domain_size, method='srt',stencil='D3Q19',
            relaxation_rate=self.omega, periodicity=(True, False, False),
            optimization={'target':self.cfg.get("macroscopic-solver", "optimization-target") , 
            'gpu_indexing':'line', 
            'double_precision':self.cfg.get("macroscopic-solver", "double-precision")})
        def obstacle(x, y, z):
            return (x > 0.45*self.cpm) & (x < 0.55*self.cpm) & (y > 0.15*self.cpm) & (y < 0.25*self.cpm)
        wall = NoSlip()
        self.scen.boundary_handling.set_boundary(wall, mask_callback=obstacle)

        def velocity_info_callback(boundary_data, **_):
            boundary_data['vel_1'] = 0
            boundary_data['vel_2'] = 0
            u_max = self.scaling_result.lattice_velocity
            y, z = boundary_data.link_positions(1), boundary_data.link_positions(2)
            H = self.domain_size[1]
            boundary_data['vel_0'] = 16 * u_max * y * z * (H-y) * (H-z) / (H*H*H*H)
        inflow = UBB(velocity_info_callback, dim=self.scen.method.dim)
        self.scen.boundary_handling.set_boundary(inflow, make_slice[0,:,:])

        outflow = FixedDensity(1.0)
        self.scen.boundary_handling.set_boundary(outflow, make_slice[-1, :, :])

    def advance(self, timesteps):
        self.scen.run(timesteps)
        self.timesteps_finished = self.timesteps_finished + timesteps
        #self.scen.write_vtk()
        self.compute_drag_lift()

    def compute_drag_lift(self):
        #scaling
        dx = self.scaling.dx
        dt = self.scaling_result.dt

        S = 0.1 * 0.41     # surface of one cylinder face in m^2 

        # compute cell indices around cylinder obstacle
        # TODO only tested for cpm = 312 cells per meter
        left = int(0.45*self.cpm) - 1  # X: left = fluid,   left+1 = obstacle
        right = int(0.55*self.cpm) + 1 # X: right = fluid,  right-1 = obstacle
        bottom = int(0.15*self.cpm)    # Y: bottom = fluid, bottom+1 = obstacle
        top = int(0.25*self.cpm)       # Y: top = fluid,    top-1 = obstacle

        # compute drag
        kvs = self.scen
        # compute advective drag force
        # pressure = density * c_s^2 = density / 3
        p_left = np.mean(kvs.density[left, bottom+1:top-1, :]) / 3
        p_right = np.mean(kvs.density[right, bottom+1:top-1, :]) / 3
        # convert from lattice to SI units
        p_left = p_left *dx*dx/(dt*dt)
        p_right = p_right *dx*dx/(dt*dt)
        FD = (p_left - p_right) * S
        # compute viscous drag force
        # ignore 3 cell layers, because simple wall boundary condition resets all PDFs to equilibrium => wrong results very close to boundary
        # use second order forward finite difference to get first derivative in normal direction of tangential velocity
        vel_grad_bottom = \
        -3/2 * np.mean(kvs.velocity[left+1:right-1, bottom-3, :, 0]) + \
        2    * np.mean(kvs.velocity[left+1:right-1, bottom-4, :, 0]) + \
        -1/2 * np.mean(kvs.velocity[left+1:right-1, bottom-5, :, 0]) 
        vel_grad_top = \
        -3/2 * np.mean(kvs.velocity[left+1:right-1, top+3, :, 0]) + \
        2    * np.mean(kvs.velocity[left+1:right-1, top+4, :, 0]) + \
        -1/2 * np.mean(kvs.velocity[left+1:right-1, top+5, :, 0])
        lb_log.debug("vel_grad_bottom = " + str(vel_grad_bottom))
        lb_log.debug("vel_grad_top = " + str(vel_grad_top))
        #convert to SI, integrate over surface, add to force
        FD = FD + self.vis * (vel_grad_bottom + vel_grad_top) / dt * S
        
        # compute lift
        # compute advective lift force
        p_bottom = np.mean(kvs.density[left+1:right-1, bottom, :]) / 3
        p_top = np.mean(kvs.density[left+1:right-1, top, :]) / 3
        p_bottom = p_bottom *dx*dx/(dt*dt)
        p_top = p_top *dx*dx/(dt*dt)
        FL = (p_bottom - p_top) * S
        # compute viscous lift force
        vel_grad_left = \
        -3/2 * np.mean(kvs.velocity[left-3, bottom+1:top-1, :, 1]) + \
        2    * np.mean(kvs.velocity[left-4, bottom+1:top-1, :, 1]) + \
        -1/2 * np.mean(kvs.velocity[left-5, bottom+1:top-1, :, 1])
        vel_grad_right = \
        -3/2 * np.mean(kvs.velocity[right+3, bottom+1:top-1, :, 1]) + \
        2    * np.mean(kvs.velocity[right+4, bottom+1:top-1, :, 1]) + \
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
        lb_log.info("Finished " + str(self.timesteps_finished) + " LB timesteps")
        lb_log.info("cD_max = " + str(self.cD_max))
        lb_log.info("cL_max = " + str(self.cL_max))

    # write velocity slice to png file
    def plot(self, filename):
        plt.vector_field_magnitude(self.scen.velocity[:,:,int(self.domain_size[2]//2),0:2])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)
        plt.axis('off')
        plt.savefig(filename)

# Read config file, create a KVSTest instance and run it
def main():
    cfg = ConfigParser()
    cfg.read("kvstest.ini")

    coloredlogs.install(fmt=
        '%(asctime)s.%(msecs)03d %(name)s %(levelname)s %(message)s'
    , level='DEBUG')

    log.setLevel(level=logging.INFO)
    lb_log.setLevel(level=logging.INFO)

    test = KVSTest(cfg)
    test.run()

# only if this file is executed directly ie. not imported as a module,
# then actually run the simulation
if __name__ == '__main__':
     main()