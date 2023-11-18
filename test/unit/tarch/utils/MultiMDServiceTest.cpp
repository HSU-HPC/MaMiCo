#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "tarch/utils/MultiMDService.h"
#if (TARCH_PARALLEL == TARCH_YES)
#include <mpi.h>
#endif

using namespace tarch;
using namespace utils;

class MultiMDServiceTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(MultiMDServiceTest);
#if (TARCH_PARALLEL == TARCH_YES)
  CPPUNIT_TEST(testAvgNumberOfMDSimulations);
  CPPUNIT_TEST(testNumberLocalComms);
  CPPUNIT_TEST(testLocalNumberOfMDSimulations);
#endif
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    _multiMDService_122_1 = new MultiMDService<3>(tarch::la::Vector<3, unsigned int>{1,2,2}, 1);
    _multiMDService_122_2 = new MultiMDService<3>(tarch::la::Vector<3, unsigned int>{1,2,2}, 2);
    _multiMDService_121_2 = new MultiMDService<3>(tarch::la::Vector<3, unsigned int>{1,2,1}, 2);
    _multiMDService_111_4 = new MultiMDService<3>(tarch::la::Vector<3, unsigned int>{1,1,1}, 4);
    _multiMDService_111_10 = new MultiMDService<3>(tarch::la::Vector<3, unsigned int>{1,1,1}, 10);
  }

  void tearDown() {
    delete _multiMDService_122_1;
    delete _multiMDService_122_2;
    delete _multiMDService_121_2;
    delete _multiMDService_111_4;
    delete _multiMDService_111_10;
  }

  void testAvgNumberOfMDSimulations() {
    unsigned int expectedAvgNumberOfMDSimulations = 1;
    unsigned int actualAvgNumberOfMDSimulations = _multiMDService_122_1->getAvgNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedAvgNumberOfMDSimulations, actualAvgNumberOfMDSimulations);

    expectedAvgNumberOfMDSimulations = 2;
    actualAvgNumberOfMDSimulations = _multiMDService_122_2->getAvgNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedAvgNumberOfMDSimulations, actualAvgNumberOfMDSimulations);

    expectedAvgNumberOfMDSimulations = 1;
    actualAvgNumberOfMDSimulations = _multiMDService_121_2->getAvgNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedAvgNumberOfMDSimulations, actualAvgNumberOfMDSimulations);

    expectedAvgNumberOfMDSimulations = 1;
    actualAvgNumberOfMDSimulations = _multiMDService_111_4->getAvgNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedAvgNumberOfMDSimulations, actualAvgNumberOfMDSimulations);

    expectedAvgNumberOfMDSimulations = 2;
    actualAvgNumberOfMDSimulations = _multiMDService_111_10->getAvgNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedAvgNumberOfMDSimulations, actualAvgNumberOfMDSimulations);
  }

  void testNumberLocalComms() {
    unsigned int expectedNumberOfLocalComms = 1;
    unsigned int actualNumberOfLocalComms = _multiMDService_122_1->getNumberLocalComms();
    CPPUNIT_ASSERT_EQUAL(expectedNumberOfLocalComms, actualNumberOfLocalComms);

    expectedNumberOfLocalComms = 1;
    actualNumberOfLocalComms = _multiMDService_122_2->getNumberLocalComms();
    CPPUNIT_ASSERT_EQUAL(expectedNumberOfLocalComms, actualNumberOfLocalComms);

    expectedNumberOfLocalComms = 2;
    actualNumberOfLocalComms = _multiMDService_121_2->getNumberLocalComms();
    CPPUNIT_ASSERT_EQUAL(expectedNumberOfLocalComms, actualNumberOfLocalComms);

    expectedNumberOfLocalComms = 4;
    actualNumberOfLocalComms = _multiMDService_111_4->getNumberLocalComms();
    CPPUNIT_ASSERT_EQUAL(expectedNumberOfLocalComms, actualNumberOfLocalComms);

    expectedNumberOfLocalComms = 4;
    actualNumberOfLocalComms = _multiMDService_111_10->getNumberLocalComms();
    CPPUNIT_ASSERT_EQUAL(expectedNumberOfLocalComms, actualNumberOfLocalComms);
  }

  void testLocalNumberOfMDSimulations() {
    CPPUNIT_ASSERT_EQUAL(_multiMDService_122_1->getLocalNumberOfMDSimulations(), _multiMDService_122_1->getAvgNumberOfMDSimulations());

    CPPUNIT_ASSERT_EQUAL(_multiMDService_122_2->getLocalNumberOfMDSimulations(), _multiMDService_122_2->getAvgNumberOfMDSimulations());

    CPPUNIT_ASSERT_EQUAL(_multiMDService_121_2->getLocalNumberOfMDSimulations(), _multiMDService_121_2->getAvgNumberOfMDSimulations());

    CPPUNIT_ASSERT_EQUAL(_multiMDService_111_4->getLocalNumberOfMDSimulations(), _multiMDService_111_4->getAvgNumberOfMDSimulations());

    unsigned int expectedLocalNumberOfMDSimulations = 2;
    int rank = 0;
  #if (TARCH_PARALLEL == TARCH_YES)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    if (rank == 3) {
      expectedLocalNumberOfMDSimulations = _multiMDService_111_10->getTotalNumberOfMDSimulations() - expectedLocalNumberOfMDSimulations * (_multiMDService_111_10->getNumberLocalComms()-1);
    }
    unsigned int actualLocalNumberOfMDSimulations = _multiMDService_111_10->getLocalNumberOfMDSimulations();
    CPPUNIT_ASSERT_EQUAL(expectedLocalNumberOfMDSimulations, actualLocalNumberOfMDSimulations);
  }

private:
  MultiMDService<3> *_multiMDService_122_1, *_multiMDService_122_2, *_multiMDService_121_2, *_multiMDService_111_4, *_multiMDService_111_10;
};

CPPUNIT_TEST_SUITE_REGISTRATION(MultiMDServiceTest);