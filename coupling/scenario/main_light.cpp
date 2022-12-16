#include "coupling/CouplingMDDefinitions.h"
#include "coupling/scenario/LightScenario.h"
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
  runScenario(new LightScenario());

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
}