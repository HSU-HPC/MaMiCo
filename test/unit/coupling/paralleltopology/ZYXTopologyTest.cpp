#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/paralleltopology/ZYXTopology.h"

class ZYXTopologyTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(ZYXTopologyTest);
  CPPUNIT_TEST(testGetProcessCoordinates3D);
  CPPUNIT_TEST(testGetProcessCoordinates2D);
  CPPUNIT_TEST(testGetRank3D);
  CPPUNIT_TEST(testGetRank2D);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
  }

  void tearDown() {
  }

   void testGetProcessCoordinates3D() {
    unsigned int nx = 3;
    unsigned int ny = 4;
    unsigned int nz = 10;
    tarch::la::Vector<3, unsigned int> numberOfProcesses{nx, ny, nz};
    unsigned int topologyOffset = 8;
    coupling::paralleltopology::ZYXTopology<3> topology{numberOfProcesses};
    for (unsigned int rank = topologyOffset; rank < topologyOffset + nx*ny*nz; ++rank) {
      tarch::la::Vector<3, unsigned int> processCoordinates = topology.getProcessCoordinates(rank, topologyOffset);
      unsigned int x = (rank-topologyOffset)/(nz*ny);
      unsigned int y = (rank-topologyOffset-x*nz*ny)/nz;
      unsigned int z = rank-topologyOffset-y*nz-x*nz*ny;
      tarch::la::Vector<3, unsigned int> expectedProcessCoordinates{x,y,z};
      CPPUNIT_ASSERT_EQUAL(expectedProcessCoordinates, processCoordinates);
    }
  }

  void testGetProcessCoordinates2D() {
    unsigned int nx = 4;
    unsigned int ny = 5;
    tarch::la::Vector<2, unsigned int> numberOfProcesses{nx, ny};
    unsigned int topologyOffset = 3;
    coupling::paralleltopology::ZYXTopology<2> topology{numberOfProcesses};
    for (unsigned int rank = topologyOffset; rank < topologyOffset + nx*ny; ++rank) {
      tarch::la::Vector<2, unsigned int> processCoordinates = topology.getProcessCoordinates(rank, topologyOffset);
      unsigned int x = (rank-topologyOffset)/ny;
      unsigned int y = rank-topologyOffset-x*ny;
      tarch::la::Vector<2, unsigned int> expectedProcessCoordinates{x,y};
      CPPUNIT_ASSERT_EQUAL(expectedProcessCoordinates, processCoordinates);
    }
  }

  void testGetRank3D() {
    unsigned int nx = 2;
    unsigned int ny = 3;
    unsigned int nz = 5;
    tarch::la::Vector<3, unsigned int> numberOfProcesses{nx, ny, nz};
    unsigned int topologyOffset = 10;
    coupling::paralleltopology::ZYXTopology<3> topology{numberOfProcesses};
    unsigned int expectedRank=topologyOffset;
    for (unsigned int k = 0; k < nx; ++k) {
      for (unsigned int j = 0; j < ny; ++j) {
        for (unsigned int i = 0; i < nz; ++i) {
          tarch::la::Vector<3, unsigned int> processCoordinates{k,j,i};
          unsigned int rank = topology.getRank(processCoordinates, topologyOffset);
          CPPUNIT_ASSERT_EQUAL(expectedRank, rank);
          expectedRank++;
        }
      }
    }
  }

  void testGetRank2D() {
    unsigned int nx = 2;
    unsigned int ny = 3;
    tarch::la::Vector<2, unsigned int> numberOfProcesses{nx, ny};
    unsigned int topologyOffset = 5;
    coupling::paralleltopology::ZYXTopology<2> topology{numberOfProcesses};
    unsigned int expectedRank=topologyOffset;
    for (unsigned int j = 0; j < nx; ++j) {
      for (unsigned int i = 0; i < ny; ++i) {
        tarch::la::Vector<2, unsigned int> processCoordinates{j,i};
        unsigned int rank = topology.getRank(processCoordinates, topologyOffset);
        CPPUNIT_ASSERT_EQUAL(expectedRank, rank);
        expectedRank++;
      }
    }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(ZYXTopologyTest);