// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MAIN_CPP_
#define _MAIN_CPP_

#include "coupling/CouplingMDDefinitions.h"
#if (BUILD_WITH_PRECICE)
#include "coupling/scenario/precice/CouetteScenario.h"
#include "coupling/scenario/precice/EvaporationScenario.h"
#endif
#include "coupling/scenario/CouetteScenario.h"
#include <cstdlib>
#include <iostream>
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

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif

  // run scenarios
  Scenario* scenario;
#if (BUILD_WITH_PRECICE)
  if (argc > 1) {  
    std::string scenarioName(argv[1]);
    if (scenarioName == "couette") {
      scenario = new precice_scenario::CouetteScenario();
    } else if (scenarioName == "evaporation") {
      scenario = new precice_scenario::EvaporationScenario();
    } else { // default to couette scenario
      scenario = new CouetteScenario();    
    }
  }
#else
  scenario = new CouetteScenario();
#endif
  runScenario(scenario);

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
}

#endif // _MAIN_CPP_
