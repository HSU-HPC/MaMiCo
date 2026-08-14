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
    MainExecSpace mainExecSpace;
    std::cout << "Kokkos using execution space \"" << mainExecSpace.name() << "\" with memory space \"" << MainExecSpace::memory_space::name() << "\""
              << std::endl;
    mainExecSpace.print_configuration(std::cout);

    std::cout << "Available concurrency: " << mainExecSpace.concurrency() << std::endl;

    SimpleMDBenchSize benchsize;
    // run tests
    if(argc == 1){
      runTest(new CellIdxIterBench());
      std::cout << std::endl << "==================== ==================== ====================" << std::endl << std::endl;
      benchsize = SimpleMDBenchSize::MD60;
    } else if(argc == 2 && !strcmp(argv[1],"MD60")) benchsize = SimpleMDBenchSize::MD60;
      else if(argc == 2 && !strcmp(argv[1],"MD120")) benchsize = SimpleMDBenchSize::MD120;
      else if(argc == 2 && !strcmp(argv[1],"MD240")) benchsize = SimpleMDBenchSize::MD240;
      else if(argc == 2 && !strcmp(argv[1],"MD480")) benchsize = SimpleMDBenchSize::MD480;
    else {
      std::cout << "ERROR unknown parameter" << std::endl;
      Kokkos::finalize();
      exit(EXIT_FAILURE);
    }
    runTest(new SimpleMDBench(benchsize));
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
