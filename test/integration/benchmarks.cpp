#include "coupling/CouplingMDDefinitions.h"
#include "test/integration/Test.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include "test/integration/CellIdxIterBench.h"
#include "test/integration/SimpleMDBench.h"

#include <Kokkos_Core.hpp>

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
  {
    Kokkos::ScopeGuard kokkos(argc, argv);
    std::cout << "Kokkos using execution space \"" << MainExecSpace::name() << "\" with memory space \"" << MainExecSpace::memory_space::name() << "\""
              << std::endl;
    MainExecSpace().print_configuration(std::cout);

    std::cout << "Available concurrency: " << MainExecSpace::concurrency() << std::endl;

    // run tests
    runTest(new CellIdxIterBench());
    std::cout << std::endl << "==================== ==================== ====================" << std::endl << std::endl;
    runTest(new SimpleMDBench());
  }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif

  return 0;
};

/* Sample CellIdxIterBench Output:

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