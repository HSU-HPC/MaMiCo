#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "tarch/la/ScalarOperations.h"

class ScalarOperationsTest : public CppUnit::TestCase
{

  CPPUNIT_TEST_SUITE(ScalarOperationsTest);
  CPPUNIT_TEST_SUITE_END();

public: 
  ScalarOperationsTest() : CppUnit::TestCase("scalar operations test") {};

  void runTest() override {
    CPPUNIT_ASSERT(tarch::la::equals<double>(2.0, 2.0, 1e-6));
    CPPUNIT_ASSERT(!tarch::la::equals<double>(2.0, 3.0, 1e-6));
  };

};

 CPPUNIT_TEST_SUITE_REGISTRATION(ScalarOperationsTest);