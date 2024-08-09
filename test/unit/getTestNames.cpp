#include "coupling/CouplingMDDefinitions.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif
#include <cppunit/Test.h>
#include <cppunit/TestSuite.h>
#include <cppunit/extensions/TestFactoryRegistry.h>

/***
 *
 * Prints testPath of all tests to stdout, for use by cmake
 *
 * @author Piet
 * Nov 2023
 *
 */
int main(int argc, char** argv) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Init(&argc, &argv);
#endif

  CppUnit::TestSuite* all = dynamic_cast<CppUnit::TestSuite*>(CppUnit::TestFactoryRegistry::getRegistry().makeTest());
  auto suites = all->getTests();
  for (auto s : suites) {
    auto suite = dynamic_cast<CppUnit::TestSuite*>(s);
    auto tests = suite->getTests();
    for (auto test : tests)
      std::cout << test->getName() << " ";
  }

  std::cout << std::endl;

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Finalize();
#endif
  return 0;
}