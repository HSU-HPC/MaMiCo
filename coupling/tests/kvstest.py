#!/usr/bin/env python3

# This file is part of the Mamico project
# Author: Piet Jarmatz
# July 2020
# Helmut Schmidt University, Hamburg. Chair for High Performance Computing
# BSD license, see the copyright notice in Mamico's main folder

import sys
import logging
import math
import mamico.tarch.configuration
from mamico.coupling.services import MultiMDCellService
from mamico.coupling.solvers import CouetteSolverInterface
from mamico.coupling.interface import MamicoInterfaceProvider
import mamico.tarch.utils
import numpy as np
from configparser import ConfigParser
#import matplotlib.animation as animation
#from matplotlib import pyplot as matplt
from mamico.coupling import getMDSimulation
from mamico.coupling import getMDSolverInterface

log = logging.getLogger('KVSTest')

# Versatile configurable MPI parallel Kármán vortex street flow test for noise-filtered multi-instance Nie coupling.
# Features:
# -> one-way or two-way coupling (TODO)
# -> using lbmpy for Lattice Bolzmann flow simulation with CUDA on GPU
# -> using (multi-instance) SimpleMD or synthetic MD data (with Gaussian noise) (TODO)
# -> using filtering subsystem with configurable coupling data analysis or noise filter sequences (TODO)
# -> separate runtime measurements for coupled simulation components (TODO)
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
        self.shutdown()

    def runOneCouplingCycle(self, cycle):
        self.advanceMacro(cycle)
        self.advanceMicro(cycle)
        self.twoWayCoupling(cycle)
        if self.rank==0 and cycle%20==0:
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
        self.macroscopicSolver = LBSolver(self.cfg)

        self.multiMDService = mamico.tarch.utils.MultiMDService(numberProcesses = self.simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(),
            totalNumberMDSimulations = self.cfg.getint("microscopic-solver","number-md-simulations"))

        self.localMDInstances = self.multiMDService.getLocalNumberOfMDSimulations()
        self.simpleMD = [getMDSimulation(self.simpleMDConfig,self.mamicoConfig, 
            self.multiMDService.getLocalCommunicator()) for i in range(self.localMDInstances)]

        for i in range(self.localMDInstances):
            self.simpleMD[i].init(self.multiMDService, self.multiMDService.getGlobalNumberOfLocalMDSimulation(i))
        
        equSteps = self.cfg.getint("microscopic-solver", "equilibration-steps")
        for i in range(self.localMDInstances):
            self.simpleMD[i].switchOffCoupling()
            self.simpleMD[i].simulateTimesteps(equSteps,0)
            self.simpleMD[i].switchOnCoupling()
        self.mdStepCounter = equSteps
    
        self.mdSolverInterface = [getMDSolverInterface(self.simpleMDConfig, self.mamicoConfig,
            self.simpleMD[i])for i in range(self.localMDInstances)]

        self.macroscopicSolverInterface = CouetteSolverInterface(self.getGlobalNumberMacroscopicCells(),
            self.mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap())
        MamicoInterfaceProvider.setMacroscopicSolverInterface(self.macroscopicSolverInterface)
    
        self.multiMDCellService = MultiMDCellService(self.mdSolverInterface, self.macroscopicSolverInterface, 
            self.simpleMDConfig, self.rank, self.cfg.getint("microscopic-solver","number-md-simulations"),
            self.mamicoConfig,self.multiMDService)

        for i in range(self.localMDInstances):
            self.simpleMD[i].setMacroscopicCellService(self.multiMDCellService.getMacroscopicCellService(i))
            self.multiMDCellService.getMacroscopicCellService(i).computeAndStoreTemperature(
                self.cfg.getfloat("microscopic-solver", "temperature"))
      
        #todo
        #allocateSendBuffer(self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface)
        #allocateRecvBuffer(self.multiMDCellService.getMacroscopicCellService(0).getIndexConversion(),*couetteSolverInterface)
    
        if self.rank==0:
            log.info("Finished initSolvers") # after ? ms

    def shutdown(self):
        pass #todo

    def advanceMacro(self,cycle):
        self.macroscopicSolver.advance(dt) #todo get dt
        pass

    def advanceMicro(self, cycle):
        for i in range(self.localMDInstances):
            self.simpleMD[i].simulateTimesteps(dt,cycle)  #todo get dt here

    def twoWayCoupling(self, cycle):
        pass #todo

    def __del__(self):
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
        self.scaling_result = self.scaling.diffusive_scaling(self.omega)
        lb_log.info(self.scaling_result)
        lb_log.info("dx=" + str(self.scaling.dx))

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

    def advance(self, dt):
        timesteps = dt / self.scaling_result.dt
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
        left = int(0.45*cpm) - 1  # X: left = fluid,   left+1 = obstacle
        right = int(0.55*cpm) + 1 # X: right = fluid,  right-1 = obstacle
        bottom = int(0.15*cpm)    # Y: bottom = fluid, bottom+1 = obstacle
        top = int(0.25*cpm)       # Y: top = fluid,    top-1 = obstacle

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
        FD = FD + vis * (vel_grad_bottom + vel_grad_top) / dt * S
        
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
        FL = FL + vis * (vel_grad_left + vel_grad_right) / dt * S

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

# Read config file, create a KVSTest instance and run it
def main():
    cfg = ConfigParser()
    cfg.read("kvstest.cfg")

    logging.basicConfig(level=logging.INFO)

    test = KVSTest(cfg)
    test.run()

# only if this file is executed directly ie. not imported as a module,
# then actually run the simulation
if __name__ == '__main__':
     main()
