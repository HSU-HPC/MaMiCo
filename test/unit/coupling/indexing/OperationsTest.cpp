#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/indexing/IndexingService.h"

#include "coupling/indexing/Operations.h"

class OperationsTest : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(OperationsTest);
	CPPUNIT_TEST_EXCEPTION(testConvertToScalarOutOfDomain, std::runtime_error);
	CPPUNIT_TEST( testConvertToScalarInDomain );
	CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {

  }


	void testConvertToScalarOutOfDomain(){
		CPPUNIT_ASSERT( true );
		throw std::runtime_error("blabla");
	}

	void testConvertToScalarInDomain(){
		CPPUNIT_ASSERT( true );
	}

private:


};

CPPUNIT_TEST_SUITE_REGISTRATION(OperationsTest);