// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MAIN_CPP_
#define _MOLECULARDYNAMICS_MAIN_CPP_

#include "simplemd/MolecularDynamicsDefinitions.h"
#if (MD_PARALLEL == MD_YES)
#include <mpi.h>
#endif
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

int main(int argc, char* argv[]) {
// initialise parallel environment
#if (MD_PARALLEL == MD_YES)
  MPI_Init(&argc, &argv);
#endif
  if (argc != 2) {
    std::cout << "Wrong number of arguments! Please call program by "
                 "./moleculardynamics inputfile.xml!"
              << std::endl;
    return -1;
  }

  // parse configuration
  simplemd::configurations::MolecularDynamicsConfiguration configuration;
  const std::string filename(argv[1]);
  tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename, "molecular-dynamics",
                                                                                                                         configuration);
  if (!configuration.isValid()) {
    std::cout << "Unvalid configuration!" << std::endl;
    return -1;
  }

  // initialise MD simulation
  simplemd::MolecularDynamicsSimulation simulation(configuration);
  simulation.initServices();

  // solve MD simulation
  for (unsigned int t = 0; t < configuration.getSimulationConfiguration().getNumberOfTimesteps(); t++) {
    simulation.simulateOneTimestep(t);
  }

  simulation.shutdownServices();

// shutdown parallel environment
#if (MD_PARALLEL == MD_YES)
  MPI_Finalize();
#endif
  return 0;
}

#endif
