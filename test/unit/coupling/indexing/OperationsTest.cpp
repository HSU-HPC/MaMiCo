#include "coupling/CouplingMDDefinitions.h"
#include "coupling/indexing/IndexingService.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace coupling::indexing;

class OperationsTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(OperationsTest);
  CPPUNIT_TEST(testConvertToScalarOutOfDomain);
  CPPUNIT_TEST(testConvertToScalarIsScalar);
  CPPUNIT_TEST(testConvertToScalarInDomain);
  CPPUNIT_TEST(testConvertTwice);
  CPPUNIT_TEST(testConvertToVectorIsVector);
  CPPUNIT_TEST_SUITE_END();

public:
  using Trait = IndexTrait;
  using Global_IdxT = CellIndex<3, Trait::vector>;
  using Global_IdxT_s = CellIndex<3>;
  using LocalMD2M_IdxT = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro, Trait::noGhost>;
  using LocalMD2M_IdxT_s = CellIndex<3, Trait::local, Trait::md2macro, Trait::noGhost>;

  void setUp() {
    _rank = 0;
    _size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

    tarch::la::Vector<3, unsigned int> numberProcesses{(unsigned int)_size, 1, 1};
    if (_size == 4)
      numberProcesses = {2, 2, 1};

    IndexingService<3>::getInstance().initWithCells({12}, numberProcesses, {1}, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
  }

  void tearDown() { IndexingService<3>::getInstance().finalize(); }

  void testConvertToScalarOutOfDomain() {
    // Contains indices of cells outside the domain
    std::vector<Global_IdxT> global_corners{{{-1, -1, -1}}, {{14, -1, -1}}, {{-1, 14, -1}}, {{14, 14, -1}},
                                            {{-1, -1, 14}}, {{14, -1, 14}}, {{-1, 14, 14}}, {{14, 14, 14}}};
    for (const auto& corner : global_corners)
      CPPUNIT_ASSERT_THROW(convertToScalar(corner), std::runtime_error);
    std::vector<Global_IdxT> md2m_corners{{{3, 3, 3}}, {{10, 3, 3}}, {{3, 10, 3}}, {{10, 10, 3}}, {{3, 3, 10}}, {{10, 3, 10}}, {{3, 10, 10}}, {{10, 10, 10}}};
    for (const auto& corner : md2m_corners) {
      LocalMD2M_IdxT local_index{corner};
      CPPUNIT_ASSERT_THROW(convertToScalar(local_index), std::runtime_error);
    }
  }

  void testConvertToScalarIsScalar() {
    using Idx_T = CellIndex<3, Trait::local, Trait::md2macro>;
    for (const auto& idx : Idx_T{})
      CPPUNIT_ASSERT_EQUAL(idx.get(), convertToScalar(idx));
  }

  void testConvertToScalarInDomain() {
    std::vector<Global_IdxT> global_corners{{{0, 0, 0}}, {{13, 0, 0}}, {{0, 13, 0}}, {{13, 13, 0}}, {{0, 0, 13}}, {{13, 0, 13}}, {{0, 13, 13}}, {{13, 13, 13}}};
    for (const auto& corner : global_corners)
      CPPUNIT_ASSERT_NO_THROW(convertToScalar(corner));
  }

  void testConvertTwice() {
    for (const auto& idx : Global_IdxT{})
      CPPUNIT_ASSERT_EQUAL(convertToVector(Global_IdxT_s{convertToScalar(idx)}), idx.get());

    for (const auto& idx : Global_IdxT_s{})
      CPPUNIT_ASSERT_EQUAL(convertToScalar(Global_IdxT{convertToVector(idx)}), idx.get());

    for (const auto& idx : LocalMD2M_IdxT{})
      CPPUNIT_ASSERT_EQUAL(convertToVector(LocalMD2M_IdxT_s{convertToScalar(idx)}), idx.get());

    for (const auto& idx : LocalMD2M_IdxT_s{})
      CPPUNIT_ASSERT_EQUAL(convertToScalar(LocalMD2M_IdxT{convertToVector(idx)}), idx.get());
  }

  void testConvertToVectorIsVector() {
    using Idx_T = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro>;
    for (const auto& idx : Idx_T{})
      CPPUNIT_ASSERT_EQUAL(idx.get(), convertToVector(idx));
  }

private:
  int _size, _rank;
};

CPPUNIT_TEST_SUITE_REGISTRATION(OperationsTest);