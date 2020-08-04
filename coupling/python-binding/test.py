#!/usr/bin/env python3

import sys
import mamico.tarch.utils
import mamico.tarch.configuration

from mamico.coupling import getMDSimulation

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

multiMDService = mamico.tarch.utils.MultiMDService(totalNumberMDSimulations = 12, 
	numberProcesses = simpleMDConfig.getMPIConfiguration().getNumberOfProcesses())

localMDInstances = multiMDService.getLocalNumberOfMDSimulations()
print("localMDInstances = " + str(localMDInstances))

print(multiMDService.getLocalCommunicator())

simpleMD = [getMDSimulation(simpleMDConfig,mamicoConfig, 
            multiMDService.getLocalCommunicator()) for i in range(localMDInstances)]

for i in range(localMDInstances):
    simpleMD[i].init(multiMDService, multiMDService.getGlobalNumberOfLocalMDSimulation(i))

mamico.tarch.utils.finalizeMPI()
