#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/molecule-mappings/VelocityStoermerVerletMapping.h"
#include "coupling/CouplingMDDefinitions.h"

#include <sstream>

class VelocityStoermerVerletMappingTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(VelocityStoermerVerletMappingTest);
  CPPUNIT_TEST(testHandleMolecule);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    
  }
  void tearDown() {}
  void testHandleMolecule() {
	
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(VelocityStoermerVerletMappingTest);