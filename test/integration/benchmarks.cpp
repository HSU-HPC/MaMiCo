#include "coupling/CouplingMDDefinitions.h"
#include "test/integration/Test.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include "test/integration/CellIdxIterBench.h"

void runTest(Test* test) {
  if (test == NULL) {
    std::cout << "ERROR executeTest: test==NULL!" << std::endl;
    exit(EXIT_FAILURE);
  }
  test->run();
  delete test;
}

int main(int argc, char* argv[]) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif

  // run tests
  runTest(new CellIdxIterBench());

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
};

/* Sampe Output:

Run CellIdxIterBench...
Number cells in test domain: 74088
lowerBoundary = 4 , 4 , 4
upperBoundary = 45 , 45 , 45

Scalar benchmark -------------
Useless result: 27444788280000
Raw loop: 191ms
Useless result: 27444788280000
Index range iterator: 119ms

Vector benchmark -------------
Useless result: 15188040000 , 15188040000 , 15188040000
Raw loop: 147ms
Useless result: 15188040000 , 15188040000 , 15188040000
Index range iterator: 152ms
Shut down CellIdxIterBench

*/