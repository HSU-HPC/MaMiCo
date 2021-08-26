// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MAIN_CPP_
#define _MOLECULARDYNAMICS_MAIN_CPP_

#include <vector>
#include "simplemd/MolecularDynamicsDefinitions.h"
#if (MD_PARALLEL==MD_YES)
#include <mpi.h>
#endif
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"

int main(int argc, char *argv[]){
  // initialise parallel environment
  #if (MD_PARALLEL==MD_YES)
  MPI_Init(&argc, &argv);
  #endif
  if (argc != 3){
    std::cout << "Wrong number of arguments! Please call program by ./moleculardynamics inputfile.xml numberMDSimulations!" << std::endl;
    return -1;
  }
  // get total number of MD simulations from commandline
  const unsigned int totalNumberMDSimulations= (unsigned int) atoi(argv[2]);
  if (totalNumberMDSimulations<=0){std::cout << "ERROR! totalNumberMDSimulations <=0!" << std::endl; exit(EXIT_FAILURE);}

  // parse configuration
  simplemd::configurations::MolecularDynamicsConfiguration configuration;
  const std::string filename(argv[1]);
  tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename,"molecular-dynamics", configuration);
  if (!configuration.isValid()){
    std::cout << "Unvalid configuration!" << std::endl;
    return -1;
  }

  // initialise multiple MD simulations
  const tarch::utils::MultiMDService<MD_DIM> multiMDService(configuration.getMPIConfiguration().getNumberOfProcesses(),totalNumberMDSimulations);
  const unsigned int localNumberMDSimulations = multiMDService.getLocalNumberOfMDSimulations();
  std::vector<simplemd::MolecularDynamicsSimulation*> simulations;
  for (unsigned int i = 0; i < localNumberMDSimulations; i++){
    simulations.push_back(new simplemd::MolecularDynamicsSimulation(configuration));
    if (simulations[i]==NULL){std::cout << "ERROR main_multi.cpp: simulations[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);}
    simulations[i]->initServices(multiMDService,multiMDService.getGlobalNumberOfLocalMDSimulation(i));
  }
  // solve MD simulations
  for (unsigned int i = 0; i < localNumberMDSimulations; i++){
    for (unsigned int t = 0; t < configuration.getSimulationConfiguration().getNumberOfTimesteps(); t++){
      simulations[i]->simulateOneTimestep(t);
    }
  }

  // shutdown MD simulations
  for (unsigned int i = 0; i < localNumberMDSimulations; i++){
    simulations[i]->shutdownServices();
    delete simulations[i]; simulations[i]=NULL;
  }
  simulations.clear();

  // shutdown parallel environment
  #if (MD_PARALLEL==MD_YES)
  MPI_Finalize();
  #endif
  return 0;
}

#endif
