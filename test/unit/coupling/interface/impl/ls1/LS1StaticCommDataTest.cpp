#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"

#include <array>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1StaticCommDataTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LS1StaticCommDataTest);
  CPPUNIT_TEST(testSerialFunctions);

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  CPPUNIT_TEST(testParallelFunctions);
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testSerialFunctions() {
    std::string filename = "test";
    coupling::interface::LS1StaticCommData::getInstance().setConfigFilename(filename);
    CPPUNIT_ASSERT(coupling::interface::LS1StaticCommData::getInstance().getConfigFilename() == filename);

    std::array<double, 3> offsets({5, 5, 5});
    for (int i = 0; i < 3; i++)
      coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(i, offsets[i]);
    for (int i = 0; i < 3; i++)
      CPPUNIT_ASSERT(coupling::interface::LS1StaticCommData::getInstance().getBoxOffsetAtDim(i) == offsets[i]);
  }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  void testParallelFunctions() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank, rank, &comm);
    coupling::interface::LS1StaticCommData::getInstance().setLocalCommunicator(comm);
    CPPUNIT_ASSERT(coupling::interface::LS1StaticCommData::getInstance().getLocalCommunicator() == comm);

    std::array<int, 3> gridDec({2, 1, 4});
    for (int i = 0; i < 3; i++)
      coupling::interface::LS1StaticCommData::getInstance().setDomainGridDecompAtDim(i, gridDec[i]);
    for (int i = 0; i < 3; i++)
      CPPUNIT_ASSERT(coupling::interface::LS1StaticCommData::getInstance().getDomainGridDecompAtDim(i) == gridDec[i]);
    CPPUNIT_ASSERT(coupling::interface::LS1StaticCommData::getInstance().getDomainGridDecomp()[0] == gridDec[0] &&
                   coupling::interface::LS1StaticCommData::getInstance().getDomainGridDecomp()[1] == gridDec[1] &&
                   coupling::interface::LS1StaticCommData::getInstance().getDomainGridDecomp()[2] == gridDec[2]);
    MPI_Comm_free(&comm);
  }
#endif
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1StaticCommDataTest);