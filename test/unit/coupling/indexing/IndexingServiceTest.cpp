// For testing private methods of IndexingService
#include "coupling/indexing/IndexingService.h"
#include "coupling/CouplingMDDefinitions.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <functional>

using namespace coupling::indexing;

class IndexingServiceTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(IndexingServiceTest);
  CPPUNIT_TEST(testAllBoundaries);
  CPPUNIT_TEST(testIndexingServiceMustBeInitialized);
  CPPUNIT_TEST(testEdgeCases);
  CPPUNIT_TEST_SUITE_END();

public:
  using Trait = IndexTrait;

  using T0 = CellIndex<3>;
  using T1 = CellIndex<3, Trait::vector>;
  using T2 = CellIndex<3, Trait::local>;
  using T3 = CellIndex<3, Trait::vector, Trait::local>;
  using T4 = CellIndex<3, Trait::md2macro>;
  using T5 = CellIndex<3, Trait::vector, Trait::md2macro>;
  using T6 = CellIndex<3, Trait::local, Trait::md2macro>;
  using T7 = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro>;
  using T8 = CellIndex<3, Trait::noGhost>;
  using T9 = CellIndex<3, Trait::vector, Trait::noGhost>;
  using T10 = CellIndex<3, Trait::local, Trait::noGhost>;
  using T11 = CellIndex<3, Trait::vector, Trait::local, Trait::noGhost>;
  using T12 = CellIndex<3, Trait::md2macro, Trait::noGhost>;
  using T13 = CellIndex<3, Trait::vector, Trait::md2macro, Trait::noGhost>;
  using T14 = CellIndex<3, Trait::local, Trait::md2macro, Trait::noGhost>;
  using T15 = CellIndex<3, Trait::vector, Trait::local, Trait::md2macro, Trait::noGhost>;

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

  void testAllBoundaries() {
    // "golden master test" for standard domain with 12x12x12 cells (plus ghost) on 4 MPI ranks

    using V = tarch::la::Vector<3, int>;
    using Vu = tarch::la::Vector<3, unsigned int>;

    CPPUNIT_ASSERT_EQUAL(T0::lowerBoundary.get(), V(0, 0, 0));
    CPPUNIT_ASSERT_EQUAL(T0::upperBoundary.get(), V(13, 13, 13));
    CPPUNIT_ASSERT_EQUAL(T0::numberCellsInDomain, Vu(14, 14, 14));
    CPPUNIT_ASSERT_EQUAL(T0::linearNumberCellsInDomain, 2744u);

    CPPUNIT_ASSERT_EQUAL(T1::lowerBoundary.get(), V(0, 0, 0));
    CPPUNIT_ASSERT_EQUAL(T1::upperBoundary.get(), V(13, 13, 13));
    CPPUNIT_ASSERT_EQUAL(T1::numberCellsInDomain, Vu(14, 14, 14));
    CPPUNIT_ASSERT_EQUAL(T1::linearNumberCellsInDomain, 2744u);

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T2::lowerBoundary.get(), V(0, 0, 0));
        CPPUNIT_ASSERT_EQUAL(T2::upperBoundary.get(), V(7, 7, 13));
        CPPUNIT_ASSERT_EQUAL(T2::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T2::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T2::lowerBoundary.get(), V(6, 0, 0));
        CPPUNIT_ASSERT_EQUAL(T2::upperBoundary.get(), V(13, 7, 13));
        CPPUNIT_ASSERT_EQUAL(T2::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T2::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T2::lowerBoundary.get(), V(0, 6, 0));
        CPPUNIT_ASSERT_EQUAL(T2::upperBoundary.get(), V(7, 13, 13));
        CPPUNIT_ASSERT_EQUAL(T2::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T2::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T2::lowerBoundary.get(), V(6, 6, 0));
        CPPUNIT_ASSERT_EQUAL(T2::upperBoundary.get(), V(13, 13, 13));
        CPPUNIT_ASSERT_EQUAL(T2::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T2::linearNumberCellsInDomain, 896u);
      }
    }

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T3::lowerBoundary.get(), V(0, 0, 0));
        CPPUNIT_ASSERT_EQUAL(T3::upperBoundary.get(), V(7, 7, 13));
        CPPUNIT_ASSERT_EQUAL(T3::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T3::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T3::lowerBoundary.get(), V(6, 0, 0));
        CPPUNIT_ASSERT_EQUAL(T3::upperBoundary.get(), V(13, 7, 13));
        CPPUNIT_ASSERT_EQUAL(T3::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T3::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T3::lowerBoundary.get(), V(0, 6, 0));
        CPPUNIT_ASSERT_EQUAL(T3::upperBoundary.get(), V(7, 13, 13));
        CPPUNIT_ASSERT_EQUAL(T3::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T3::linearNumberCellsInDomain, 896u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T3::lowerBoundary.get(), V(6, 6, 0));
        CPPUNIT_ASSERT_EQUAL(T3::upperBoundary.get(), V(13, 13, 13));
        CPPUNIT_ASSERT_EQUAL(T3::numberCellsInDomain, Vu(8, 8, 14));
        CPPUNIT_ASSERT_EQUAL(T3::linearNumberCellsInDomain, 896u);
      }
    }

    CPPUNIT_ASSERT_EQUAL(T4::lowerBoundary.get(), V(3, 3, 3));
    CPPUNIT_ASSERT_EQUAL(T4::upperBoundary.get(), V(10, 10, 10));
    CPPUNIT_ASSERT_EQUAL(T4::numberCellsInDomain, Vu(8, 8, 8));
    CPPUNIT_ASSERT_EQUAL(T4::linearNumberCellsInDomain, 512u);

    CPPUNIT_ASSERT_EQUAL(T5::lowerBoundary.get(), V(3, 3, 3));
    CPPUNIT_ASSERT_EQUAL(T5::upperBoundary.get(), V(10, 10, 10));
    CPPUNIT_ASSERT_EQUAL(T5::numberCellsInDomain, Vu(8, 8, 8));
    CPPUNIT_ASSERT_EQUAL(T5::linearNumberCellsInDomain, 512u);

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T6::lowerBoundary.get(), V(3, 3, 3));
        CPPUNIT_ASSERT_EQUAL(T6::upperBoundary.get(), V(7, 7, 10));
        CPPUNIT_ASSERT_EQUAL(T6::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T6::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T6::lowerBoundary.get(), V(6, 3, 3));
        CPPUNIT_ASSERT_EQUAL(T6::upperBoundary.get(), V(10, 7, 10));
        CPPUNIT_ASSERT_EQUAL(T6::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T6::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T6::lowerBoundary.get(), V(3, 6, 3));
        CPPUNIT_ASSERT_EQUAL(T6::upperBoundary.get(), V(7, 10, 10));
        CPPUNIT_ASSERT_EQUAL(T6::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T6::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T6::lowerBoundary.get(), V(6, 6, 3));
        CPPUNIT_ASSERT_EQUAL(T6::upperBoundary.get(), V(10, 10, 10));
        CPPUNIT_ASSERT_EQUAL(T6::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T6::linearNumberCellsInDomain, 200u);
      }
    }

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T7::lowerBoundary.get(), V(3, 3, 3));
        CPPUNIT_ASSERT_EQUAL(T7::upperBoundary.get(), V(7, 7, 10));
        CPPUNIT_ASSERT_EQUAL(T7::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T7::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T7::lowerBoundary.get(), V(6, 3, 3));
        CPPUNIT_ASSERT_EQUAL(T7::upperBoundary.get(), V(10, 7, 10));
        CPPUNIT_ASSERT_EQUAL(T7::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T7::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T7::lowerBoundary.get(), V(3, 6, 3));
        CPPUNIT_ASSERT_EQUAL(T7::upperBoundary.get(), V(7, 10, 10));
        CPPUNIT_ASSERT_EQUAL(T7::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T7::linearNumberCellsInDomain, 200u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T7::lowerBoundary.get(), V(6, 6, 3));
        CPPUNIT_ASSERT_EQUAL(T7::upperBoundary.get(), V(10, 10, 10));
        CPPUNIT_ASSERT_EQUAL(T7::numberCellsInDomain, Vu(5, 5, 8));
        CPPUNIT_ASSERT_EQUAL(T7::linearNumberCellsInDomain, 200u);
      }
    }

    CPPUNIT_ASSERT_EQUAL(T8::lowerBoundary.get(), V(1, 1, 1));
    CPPUNIT_ASSERT_EQUAL(T8::upperBoundary.get(), V(12, 12, 12));
    CPPUNIT_ASSERT_EQUAL(T8::numberCellsInDomain, Vu(12, 12, 12));
    CPPUNIT_ASSERT_EQUAL(T8::linearNumberCellsInDomain, 1728u);

    CPPUNIT_ASSERT_EQUAL(T9::lowerBoundary.get(), V(1, 1, 1));
    CPPUNIT_ASSERT_EQUAL(T9::upperBoundary.get(), V(12, 12, 12));
    CPPUNIT_ASSERT_EQUAL(T9::numberCellsInDomain, Vu(12, 12, 12));
    CPPUNIT_ASSERT_EQUAL(T9::linearNumberCellsInDomain, 1728u);

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T10::lowerBoundary.get(), V(1, 1, 1));
        CPPUNIT_ASSERT_EQUAL(T10::upperBoundary.get(), V(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T10::lowerBoundary.get(), V(7, 1, 1));
        CPPUNIT_ASSERT_EQUAL(T10::upperBoundary.get(), V(12, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T10::lowerBoundary.get(), V(1, 7, 1));
        CPPUNIT_ASSERT_EQUAL(T10::upperBoundary.get(), V(6, 12, 12));
        CPPUNIT_ASSERT_EQUAL(T10::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T10::lowerBoundary.get(), V(7, 7, 1));
        CPPUNIT_ASSERT_EQUAL(T10::upperBoundary.get(), V(12, 12, 12));
        CPPUNIT_ASSERT_EQUAL(T10::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T10::linearNumberCellsInDomain, 432u);
      }
    }

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T11::lowerBoundary.get(), V(1, 1, 1));
        CPPUNIT_ASSERT_EQUAL(T11::upperBoundary.get(), V(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T11::lowerBoundary.get(), V(7, 1, 1));
        CPPUNIT_ASSERT_EQUAL(T11::upperBoundary.get(), V(12, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T11::lowerBoundary.get(), V(1, 7, 1));
        CPPUNIT_ASSERT_EQUAL(T11::upperBoundary.get(), V(6, 12, 12));
        CPPUNIT_ASSERT_EQUAL(T11::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::linearNumberCellsInDomain, 432u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T11::lowerBoundary.get(), V(7, 7, 1));
        CPPUNIT_ASSERT_EQUAL(T11::upperBoundary.get(), V(12, 12, 12));
        CPPUNIT_ASSERT_EQUAL(T11::numberCellsInDomain, Vu(6, 6, 12));
        CPPUNIT_ASSERT_EQUAL(T11::linearNumberCellsInDomain, 432u);
      }
    }

    CPPUNIT_ASSERT_EQUAL(T12::lowerBoundary.get(), V(4, 4, 4));
    CPPUNIT_ASSERT_EQUAL(T12::upperBoundary.get(), V(9, 9, 9));
    CPPUNIT_ASSERT_EQUAL(T12::numberCellsInDomain, Vu(6, 6, 6));
    CPPUNIT_ASSERT_EQUAL(T12::linearNumberCellsInDomain, 216u);

    CPPUNIT_ASSERT_EQUAL(T13::lowerBoundary.get(), V(4, 4, 4));
    CPPUNIT_ASSERT_EQUAL(T13::upperBoundary.get(), V(9, 9, 9));
    CPPUNIT_ASSERT_EQUAL(T13::numberCellsInDomain, Vu(6, 6, 6));
    CPPUNIT_ASSERT_EQUAL(T13::linearNumberCellsInDomain, 216u);

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T14::lowerBoundary.get(), V(4, 4, 4));
        CPPUNIT_ASSERT_EQUAL(T14::upperBoundary.get(), V(6, 6, 9));
        CPPUNIT_ASSERT_EQUAL(T14::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T14::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T14::lowerBoundary.get(), V(7, 4, 4));
        CPPUNIT_ASSERT_EQUAL(T14::upperBoundary.get(), V(9, 6, 9));
        CPPUNIT_ASSERT_EQUAL(T14::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T14::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T14::lowerBoundary.get(), V(4, 7, 4));
        CPPUNIT_ASSERT_EQUAL(T14::upperBoundary.get(), V(6, 9, 9));
        CPPUNIT_ASSERT_EQUAL(T14::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T14::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T14::lowerBoundary.get(), V(7, 7, 4));
        CPPUNIT_ASSERT_EQUAL(T14::upperBoundary.get(), V(9, 9, 9));
        CPPUNIT_ASSERT_EQUAL(T14::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T14::linearNumberCellsInDomain, 54u);
      }
    }

    if (_size == 4) {
      if (_rank == 0) {
        CPPUNIT_ASSERT_EQUAL(T15::lowerBoundary.get(), V(4, 4, 4));
        CPPUNIT_ASSERT_EQUAL(T15::upperBoundary.get(), V(6, 6, 9));
        CPPUNIT_ASSERT_EQUAL(T15::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T15::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 1) {
        CPPUNIT_ASSERT_EQUAL(T15::lowerBoundary.get(), V(7, 4, 4));
        CPPUNIT_ASSERT_EQUAL(T15::upperBoundary.get(), V(9, 6, 9));
        CPPUNIT_ASSERT_EQUAL(T15::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T15::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 2) {
        CPPUNIT_ASSERT_EQUAL(T15::lowerBoundary.get(), V(4, 7, 4));
        CPPUNIT_ASSERT_EQUAL(T15::upperBoundary.get(), V(6, 9, 9));
        CPPUNIT_ASSERT_EQUAL(T15::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T15::linearNumberCellsInDomain, 54u);
      }
      if (_rank == 3) {
        CPPUNIT_ASSERT_EQUAL(T15::lowerBoundary.get(), V(7, 7, 4));
        CPPUNIT_ASSERT_EQUAL(T15::upperBoundary.get(), V(9, 9, 9));
        CPPUNIT_ASSERT_EQUAL(T15::numberCellsInDomain, Vu(3, 3, 6));
        CPPUNIT_ASSERT_EQUAL(T15::linearNumberCellsInDomain, 54u);
      }
    }
  }

  void testIndexingServiceMustBeInitialized() {
    // some operations that need the indexing system to be initialized
    std::vector<std::function<void(void)>> operations;
    tarch::la::Vector<3, int> one{1, 1, 1};
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    operations.push_back([one]() { IndexingService<3>::getInstance().getRanksForGlobalIndex(one); });
    operations.push_back([]() { IndexingService<3>::getInstance().getComm(); });
#endif
    operations.push_back([]() { IndexingService<3>::getInstance().getRank(); });
    operations.push_back([one]() { convertToScalar(T5{one}); });
    operations.push_back([]() { convertToVector(T4{1}); });
    // static_cast<void>() is used to explicity discard value (prevents unused-value warning)
    operations.push_back([one]() { static_cast<void>(T0{T5{one}}); });
    operations.push_back([]() { static_cast<void>(T5{T0{0}}); });

    // it was initialized by setup() before, so now they should be ok
    for (auto& op : operations)
      CPPUNIT_ASSERT_NO_THROW(op());

    // calls IndexingService::finalize(), now all should fail
    tearDown();
    for (auto& op : operations)
      CPPUNIT_ASSERT_THROW(op(), std::runtime_error);
  }

  void testEdgeCases() {
    auto& service = IndexingService<3>::getInstance();

    // re-initialization
    CPPUNIT_ASSERT_NO_THROW(service.initWithCells({12}, {1, 1, 1}, {1}, coupling::paralleltopology::XYZ, 3, 0u));

    // empty MD2M domain
    CPPUNIT_ASSERT_NO_THROW(service.initWithCells({6}, {1, 1, 1}, {1}, coupling::paralleltopology::XYZ, 3, 0u));

    // non divisible domain
    CPPUNIT_ASSERT_THROW(service.initWithCells({6, 6, 7}, {1, 1, 2}, {1}, coupling::paralleltopology::XYZ, 3, 0u), std::runtime_error);

    service.initWithCells({12}, {2, 1, 1}, {1}, coupling::paralleltopology::XYZ, 3, 0u);

    // out of domain
    CPPUNIT_ASSERT_THROW(service.getUniqueRankForCouplingCell({1, 14, 1}, {12}), std::runtime_error);
    CPPUNIT_ASSERT_EQUAL(service.getUniqueRankForCouplingCell({1, 1, 1}, {12}), 0u);
    CPPUNIT_ASSERT_EQUAL(service.getUniqueRankForCouplingCell({7, 1, 1}, {12}), 1u);
  }

private:
  int _size, _rank;
};

CPPUNIT_TEST_SUITE_REGISTRATION(IndexingServiceTest);