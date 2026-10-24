#ifndef _MAIN_COUETTE_CPP_
#define _MAIN_COUETTE_CPP_
#include "coupling/CouplingMDDefinitions.h"
#if (BUILD_WITH_PRECICE)
#include "coupling/scenario/precice/CouetteScenario.h"
// precice/EvaporationScenario.h omitted on purpose: it hard-depends on ls1-MarDyn
#endif
#include "coupling/scenario/CouetteScenario.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
void runScenario(Scenario* scenario) {
  if (scenario == NULL) { std::cout << "ERROR: scenario==NULL!" << std::endl; exit(EXIT_FAILURE); }
  scenario->run();
  delete scenario;
}
// usage: couette [preCICE_Couette [config.xml]]
int main(int argc, char* argv[]) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif
  Scenario* scenario = NULL;
#if (BUILD_WITH_PRECICE)
  if (argc > 1 && strcmp(argv[1], "preCICE_Couette") == 0)
    scenario = new coupling::scenario::precice::CouetteScenario(argc > 2 ? std::string(argv[2]) : std::string("couette.xml"));
  else
    scenario = new CouetteScenario();
#else
  scenario = new CouetteScenario();
#endif
  runScenario(scenario);
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif
  return 0;
}
#endif // _MAIN_COUETTE_CPP_
