#include "coupling/solvers/LBCouetteSolver.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace coupling::solvers;

/**
 *  @author Piet Jarmatz
 */
class LBCouetteSolverStateTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(LBCouetteSolverStateTest);
  CPPUNIT_TEST(testZeroInit);
  CPPUNIT_TEST(testInit);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(TestCopyAssignments);
  CPPUNIT_TEST(testOpPlus);
  CPPUNIT_TEST(testOpMinus);
  CPPUNIT_TEST(testEquality);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    double data_one[3]{1, 1, 1};
    double data_two[3]{2, 2, 2};
    _one = LBCouetteSolverState{3, data_one};
    _two = LBCouetteSolverState{3, data_two};
  }

  void tearDown() {}

  void testZeroInit() {
    LBCouetteSolverState zero(5);
    const double* data = zero.getData();
    CPPUNIT_ASSERT(data[0] == 0);
    CPPUNIT_ASSERT(data[1] == 0);
    CPPUNIT_ASSERT(data[2] == 0);
    CPPUNIT_ASSERT(data[3] == 0);
    CPPUNIT_ASSERT(data[4] == 0);
  }

  void testInit() {
    double pdf[5]{1, 2, 3, 4, 5.5};
    LBCouetteSolverState state(5, pdf);
    const double* data = state.getData();
    CPPUNIT_ASSERT(data[0] == 1);
    CPPUNIT_ASSERT(data[1] == 2);
    CPPUNIT_ASSERT(data[2] == 3);
    CPPUNIT_ASSERT(data[3] == 4);
    CPPUNIT_ASSERT(data[4] == 5.5);
  }

  void testSize() {
    double pdf[5]{1, 2, 3, 4, 5.5};
    LBCouetteSolverState state(5, pdf);
    CPPUNIT_ASSERT(state.getSizeBytes() == 40); // 5 values times 8 bytes per double
  }

  void TestCopyAssignments() {
    double pdf[3]{100, 100, 100};
    LBCouetteSolverState a(3, pdf);
    LBCouetteSolverState b = a;
    CPPUNIT_ASSERT(b.getData()[0] == 100);
    CPPUNIT_ASSERT(b.getData()[1] == 100);
    CPPUNIT_ASSERT(b.getData()[2] == 100);
    b.getData()[1] = 50;
    CPPUNIT_ASSERT(a.getData()[1] == 100);
    CPPUNIT_ASSERT(b.getData()[1] == 50);
  }

  void testOpPlus() {
    auto res1 = _one + _two;
    LBCouetteSolverState res = *dynamic_cast<LBCouetteSolverState*>(res1.get());
    CPPUNIT_ASSERT(res.getData()[0] == 3);
    CPPUNIT_ASSERT(res.getData()[1] == 3);
    CPPUNIT_ASSERT(res.getData()[2] == 3);
    auto res2 = _two + _one;
    res = *dynamic_cast<LBCouetteSolverState*>(res2.get());
    CPPUNIT_ASSERT(res.getData()[0] == 3);
    CPPUNIT_ASSERT(res.getData()[1] == 3);
    CPPUNIT_ASSERT(res.getData()[2] == 3);
  }

  void testOpMinus() {
    auto res1 = _one - _two;
    LBCouetteSolverState res = *dynamic_cast<LBCouetteSolverState*>(res1.get());
    CPPUNIT_ASSERT(res.getData()[0] == -1);
    CPPUNIT_ASSERT(res.getData()[1] == -1);
    CPPUNIT_ASSERT(res.getData()[2] == -1);
    auto res2 = _two - _one;
    res = *dynamic_cast<LBCouetteSolverState*>(res2.get());
    CPPUNIT_ASSERT(res.getData()[0] == 1);
    CPPUNIT_ASSERT(res.getData()[1] == 1);
    CPPUNIT_ASSERT(res.getData()[2] == 1);
  }

  void testEquality() {
    double data_two[3]{2, 2, 2};
    LBCouetteSolverState my_two = LBCouetteSolverState{3, data_two};
    CPPUNIT_ASSERT(my_two == _two);
    LBCouetteSolverState res1 = *dynamic_cast<LBCouetteSolverState*>((_one + _two).get());
    LBCouetteSolverState res2 = *dynamic_cast<LBCouetteSolverState*>((_two + _one).get());
    CPPUNIT_ASSERT(res1 == res2);
    CPPUNIT_ASSERT(!(my_two == res1));
    CPPUNIT_ASSERT(!(_two == res1));
    CPPUNIT_ASSERT(!(my_two == res2));
    CPPUNIT_ASSERT(!(_two == res2));
  }

private:
  LBCouetteSolverState _one{1}, _two{1};
};

CPPUNIT_TEST_SUITE_REGISTRATION(LBCouetteSolverStateTest);
