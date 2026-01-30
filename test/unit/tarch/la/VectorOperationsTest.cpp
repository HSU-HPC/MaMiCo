#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "tarch/la/Vector.h"

using namespace tarch;
using namespace la;

class VectorOperationsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(VectorOperationsTest);
  CPPUNIT_TEST(testEquality);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST(testSubstraction);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    v_1_2_3 = new Vector<3, double>(1.0, 2.0, 3.0);
    v_4_5_6 = new Vector<3, double>(4.0, 5.0, 6.0);
  }

  void tearDown() {
    delete v_1_2_3;
    delete v_4_5_6;
  }

  void testEquality() {
    CPPUNIT_ASSERT(*v_1_2_3 == *v_1_2_3);
    CPPUNIT_ASSERT(!(*v_1_2_3 == *v_4_5_6));
  }

  void testAddition() {
    Vector<3, double> v{3.0};
    CPPUNIT_ASSERT(v + *v_1_2_3 == *v_4_5_6);
  }

  void testSubstraction() {
    Vector<3, double> v{3.0};
    CPPUNIT_ASSERT(*v_4_5_6 - v == *v_1_2_3);
  }

private:
  Vector<3, double>*v_1_2_3, *v_4_5_6;
};

CPPUNIT_TEST_SUITE_REGISTRATION(VectorOperationsTest);
