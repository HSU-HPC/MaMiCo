#!/usr/bin/env python3

import sys
import math
import mamico.tarch.utils
import mamico.tarch.configuration
import mamico.coupling 
from mamico.coupling.solvers import CouetteSolverInterface

rank = mamico.tarch.utils.initMPI()
print("rank = " + str(rank))

simpleMDConfig = mamico.tarch.configuration.parseMolecularDynamicsConfiguration("../tests/build_couette/couette_simplemd.xml","molecular-dynamics")
print(simpleMDConfig)
if not simpleMDConfig.isValid():
    print("Invalid SimpleMD config!")
    sys.exit(1)

mamicoConfig = mamico.tarch.configuration.parseMaMiCoConfiguration("../tests/build_couette/couette_mamico.xml","mamico")
print(mamicoConfig)
if not mamicoConfig.isValid():
    print("Invalid MaMiCo config!")
    sys.exit(1)

multiMDService = mamico.tarch.utils.MultiMDService(totalNumberMDSimulations = 2, 
	numberProcesses = simpleMDConfig.getMPIConfiguration().getNumberOfProcesses())

localMDInstances = multiMDService.getLocalNumberOfMDSimulations()
print("localMDInstances = " + str(localMDInstances))

print(multiMDService.getLocalCommunicator())

simpleMD = [mamico.coupling.getMDSimulation(simpleMDConfig,mamicoConfig, 
            multiMDService.getLocalCommunicator()) for i in range(localMDInstances)]

for i in range(localMDInstances):
    simpleMD[i].init(multiMDService, multiMDService.getGlobalNumberOfLocalMDSimulation(i))

equSteps = 60
for i in range(localMDInstances):
    simpleMD[i].switchOffCoupling()
    simpleMD[i].simulateTimesteps(equSteps,0)
    simpleMD[i].switchOnCoupling()
mdStepCounter = equSteps

mdSolverInterface = [mamico.coupling.getMDSolverInterface(simpleMDConfig, mamicoConfig,
            simpleMD[i]) for i in range(localMDInstances)]

domainSize = simpleMDConfig.getDomainConfiguration().getGlobalDomainSize()
dx = mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()
globalNumberMacroscopicCells = [math.floor(domainSize[d]/dx[d]+0.5) for d in range(3)]

print("domainSize = " + str(domainSize))
print("dx = " + str(dx))
print("globalNumberMacroscopicCells = " + str(globalNumberMacroscopicCells))

macroscopicSolverInterface = CouetteSolverInterface(globalNumberMacroscopicCells,
    mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap())

mamico.coupling.setMacroscopicSolverInterface(macroscopicSolverInterface)

mamico.tarch.utils.finalizeMPI()
