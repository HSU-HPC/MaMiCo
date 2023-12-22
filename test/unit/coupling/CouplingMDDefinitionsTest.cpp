#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/CouplingMDDefinitions.h"

class CouplingMDDefinitionsTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(CouplingMDDefinitionsTest);
  CPPUNIT_TEST(testInitRange);
  CPPUNIT_TEST(testInitDimVector);
  CPPUNIT_TEST(testInitDivisionFactor);
  CPPUNIT_TEST(testGetVectorCellIndex);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {
  }

  void testInitRange() {
    tarch::la::Vector<3, unsigned int> vec3{2,3,4};
    tarch::la::Vector<3, unsigned int> results = coupling::initRange<3>(vec3);
    tarch::la::Vector<3, unsigned int> expected = vec3;
    CPPUNIT_ASSERT_EQUAL(expected, results);

    tarch::la::Vector<2, unsigned int> vec2{2,3};
    results = coupling::initRange<2>(vec2);
    expected = tarch::la::Vector<3, unsigned int>{vec2[0], vec2[1], 1};
    CPPUNIT_ASSERT_EQUAL(expected, results);

    tarch::la::Vector<1, unsigned int> vec1{2};
    results = coupling::initRange<1>(vec1);
    expected = tarch::la::Vector<3, unsigned int>{vec1[0], 1, 1};
    CPPUNIT_ASSERT_EQUAL(expected, results);
  }

  void testInitDimVector() {
    tarch::la::Vector<3, unsigned int> vec{2,3,4};
    tarch::la::Vector<3, unsigned int> results3 = coupling::initDimVector<3>(vec);
    tarch::la::Vector<3, unsigned int> expected3 = vec;
    CPPUNIT_ASSERT_EQUAL(expected3, results3);

    tarch::la::Vector<2, unsigned int> results2 = coupling::initDimVector<2>(vec);
    tarch::la::Vector<2, unsigned int> expected2{vec[0], vec[1]};
    CPPUNIT_ASSERT_EQUAL(expected2, results2);

    tarch::la::Vector<1, unsigned int> results1 = coupling::initDimVector<1>(vec);
    tarch::la::Vector<1, unsigned int> expected1{vec[0]};
    CPPUNIT_ASSERT_EQUAL(expected1, results1);
  }

  void testInitDivisionFactor() {
    tarch::la::Vector<3, unsigned int> vec3{2,3,4};
    tarch::la::Vector<3, unsigned int> results3 = coupling::initDivisionFactor<3>(vec3);
    tarch::la::Vector<3, unsigned int> expected3{1,2,6};
    CPPUNIT_ASSERT_EQUAL(expected3, results3);

    tarch::la::Vector<2, unsigned int> vec2{2,3};
    tarch::la::Vector<2, unsigned int> results2 = coupling::initDivisionFactor<2>(vec2);
    tarch::la::Vector<2, unsigned int> expected2{1, 2};
    CPPUNIT_ASSERT_EQUAL(expected2, results2);

    tarch::la::Vector<1, unsigned int> vec1{2};
    tarch::la::Vector<1, unsigned int> results1 = coupling::initDivisionFactor<1>(vec1);
    tarch::la::Vector<1, unsigned int> expected1{1};
    CPPUNIT_ASSERT_EQUAL(expected1, results1);
  }

  void testGetVectorCellIndex() {
    unsigned int nx = 2;
    unsigned int ny = 4;
    unsigned int nz = 5;
    tarch::la::Vector<3, unsigned int> divisionFactor3{1,nx,nx*ny};
    unsigned int cellIndex = 0;
    tarch::la::Vector<3, unsigned int> results3;
    tarch::la::Vector<3, unsigned int> expected3;
    for (unsigned int k = 0; k < nz; ++k) {
      for (unsigned int j = 0; j < ny; ++j) {
        for (unsigned int i = 0; i < nx; ++i) {
          results3 = coupling::getVectorCellIndex<3>(cellIndex, divisionFactor3);
          expected3 = tarch::la::Vector<3, unsigned int>{i,j,k};
          CPPUNIT_ASSERT_EQUAL(expected3, results3);
          cellIndex++;
        }
      }
    }

    nx = 2;
    ny = 4;
    tarch::la::Vector<2, unsigned int> divisionFactor2{1,nx};
    cellIndex = 0;
    tarch::la::Vector<2, unsigned int> results2;
    tarch::la::Vector<2, unsigned int> expected2;
    for (unsigned int j = 0; j < ny; ++j) {
      for (unsigned int i = 0; i < nx; ++i) {
        results2 = coupling::getVectorCellIndex<2>(cellIndex, divisionFactor2);
        expected2 = tarch::la::Vector<2, unsigned int>{i,j};
        CPPUNIT_ASSERT_EQUAL(expected2, results2);
        cellIndex++;
      }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(CouplingMDDefinitionsTest);