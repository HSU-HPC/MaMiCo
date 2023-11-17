#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/paralleltopology/XYZTopology.h"
#if (TARCH_PARALLEL == TARCH_YES)
#include <mpi.h>
#endif

class XYZTopologyTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(XYZTopologyTest);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {
  }


private:
};

CPPUNIT_TEST_SUITE_REGISTRATION(XYZTopologyTest);