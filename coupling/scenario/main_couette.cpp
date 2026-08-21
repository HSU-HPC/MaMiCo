// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAIN_COUETTE_CPP_
#define _MAIN_COUETTE_CPP_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/scenario/CouetteCommandLine.h"
#include "coupling/scenario/CouetteScenario.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

/** executes a newly created scenario and deletes it immediately again. */
void runScenario(Scenario* scenario) {
  if (scenario == NULL) {
    std::cout << "ERROR executeScenario: scenario==NULL!" << std::endl;
    exit(EXIT_FAILURE);
  }
  scenario->run();
  delete scenario;
}

int main(int argc, char* argv[]) {
  std::vector<std::string> arguments(argv, argv + argc);
  CouetteCommandLineOptions options;
  try {
    options = parseCouetteCommandLine(arguments);
  } catch (const std::invalid_argument& error) {
    std::cerr << "ERROR: " << error.what() << "\nUse --help for usage." << std::endl;
    return EXIT_FAILURE;
  }
  if (options.showHelp) {
    std::cout << "Usage: " << arguments.front() << " [--config <path>] [Kokkos and MPI options]\n"
              << "\n"
              << "  --config <path>  Read the Couette configuration from <path> instead of couette.xml.\n"
              << "  -h, --help       Show this help message.\n";
    return EXIT_SUCCESS;
  }

  std::vector<char*> remainingArguments;
  remainingArguments.reserve(options.remainingArguments.size() + 1);
  for (std::string& argument : options.remainingArguments)
    remainingArguments.push_back(argument.data());
  remainingArguments.push_back(nullptr);
  int remainingArgc = static_cast<int>(options.remainingArguments.size());
  char** remainingArgv = remainingArguments.data();

  int rank = 0;
  int size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&remainingArgc, &remainingArgv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  Kokkos::ScopeGuard kokkos(remainingArgc, remainingArgv);
  MainExecSpace mainExecSpace;
  if (rank == 0) {
    std::cout << std::endl;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    std::cout << "(Printing Kokkos configuration on rank " << rank << " out of " << size << ")" << std::endl;
#endif
    std::cout << "Kokkos using execution space \"" << mainExecSpace.name() << "\" with memory space \"" << MainExecSpace::memory_space::name() << "\""
              << std::endl;
    mainExecSpace.print_configuration(std::cout);

    std::cout << "Available concurrency: " << mainExecSpace.concurrency() << std::endl << std::endl;
  }

  // run scenarios
  runScenario(new CouetteScenario(options.configurationFile));

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
}

#endif // _MAIN_COUETTE_CPP_
