#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "coupling/paralleltopology/XYZTopology.h"

class XYZTopologyTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(XYZTopologyTest);
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
    unsigned int nx = 2;
    unsigned int ny = 3;
    unsigned int nz = 5;
    tarch::la::Vector<3, unsigned int> numberOfProcesses{nx, ny, nz};
    unsigned int topologyOffset = 10;
    coupling::paralleltopology::XYZTopology<3> topology{numberOfProcesses};
    for (unsigned int rank = topologyOffset; rank < topologyOffset + nx*ny*nz; ++rank) {
      tarch::la::Vector<3, unsigned int> processCoordinates = topology.getProcessCoordinates(rank, topologyOffset);
      unsigned int z = (rank-topologyOffset)/(nx*ny);
      unsigned int y = (rank-topologyOffset-z*nx*ny)/nx;
      unsigned x = rank-topologyOffset-y*nx-z*nx*ny;
      tarch::la::Vector<3, unsigned int> expectedProcessCoordinates{x,y,z};
      CPPUNIT_ASSERT_EQUAL(expectedProcessCoordinates, processCoordinates);
    }
  }

  void testGetProcessCoordinates2D() {
    unsigned int nx = 2;
    unsigned int ny = 3;
    tarch::la::Vector<2, unsigned int> numberOfProcesses{nx, ny};
    unsigned int topologyOffset = 5;
    coupling::paralleltopology::XYZTopology<2> topology{numberOfProcesses};
    for (unsigned int rank = topologyOffset; rank < topologyOffset + nx*ny; ++rank) {
      tarch::la::Vector<2, unsigned int> processCoordinates = topology.getProcessCoordinates(rank, topologyOffset);
      unsigned int y = (rank-topologyOffset)/nx;
      unsigned int x = rank-topologyOffset-y*nx;
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
    coupling::paralleltopology::XYZTopology<3> topology{numberOfProcesses};
    unsigned int expectedRank=topologyOffset;
    for (unsigned int k = 0; k < nz; ++k) {
      for (unsigned int j = 0; j < ny; ++j) {
        for (unsigned int i = 0; i < nx; ++i) {
          tarch::la::Vector<3, unsigned int> processCoordinates{i,j,k};
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
    coupling::paralleltopology::XYZTopology<2> topology{numberOfProcesses};
    unsigned int expectedRank=topologyOffset;
    for (unsigned int j = 0; j < ny; ++j) {
      for (unsigned int i = 0; i < nx; ++i) {
        tarch::la::Vector<2, unsigned int> processCoordinates{i,j};
        unsigned int rank = topology.getRank(processCoordinates, topologyOffset);
        CPPUNIT_ASSERT_EQUAL(expectedRank, rank);
        expectedRank++;
      }
    }
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(XYZTopologyTest);