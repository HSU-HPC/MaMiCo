#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

#if defined(LS1_MARDYN)
#define MY_LINKEDCELL ls1::LS1RegionWrapper
#endif

class MDSimulationFactoryTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(MDSimulationFactoryTest);
  // CPPUNIT_TEST(testSimulateTimesteps);
#if defined(LS1_MARDYN)
  //
#endif
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  // CPPUNIT_TEST(testParallelFunctions);
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
};

CPPUNIT_TEST_SUITE_REGISTRATION(MDSimulationFactoryTest);