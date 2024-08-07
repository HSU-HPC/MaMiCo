#include "coupling/datastructures/CellContainer.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <tuple>
#include <vector>

using namespace coupling::indexing;

class CellContainerTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(CellContainerTest);
  CPPUNIT_TEST(testInsertAccessSize);
  CPPUNIT_TEST(testIteration);
  CPPUNIT_TEST(testIterationLocal);
  CPPUNIT_TEST(testEquality);
  CPPUNIT_TEST(testEqualityContent);
  CPPUNIT_TEST(testInequality);
  CPPUNIT_TEST(testPreIncrement);
  CPPUNIT_TEST(testPostIncrement);
  CPPUNIT_TEST(testIndexAccess);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    using namespace coupling::datastructures;
    _rank = 0;
    _size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
    tarch::la::Vector<3, unsigned int> numberProcesses{1, 1, 1};
    if (_size == 4)
      numberProcesses = {2, 2, 1};
    IndexingService<3>::getInstance().initWithCells(8, numberProcesses, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
    int i = 0;
    _couplingCells_fullcase.reserve(I01::linearNumberCellsInDomain);
    _idxs_fullcase.reserve(I01::linearNumberCellsInDomain);
    _couplingCells_local.reserve(I03::linearNumberCellsInDomain);
    _idxs_local.reserve(I03::linearNumberCellsInDomain);
    for (auto idx : I01()) {
      auto tmp = new CouplingCell<3>;
      tmp->setMacroscopicMass(i++);
      _couplingCells_fullcase.push_back(tmp);
      _idxs_fullcase.push_back(idx);
    }
    for (auto idx : I03()) {
      _couplingCells_local.push_back(_couplingCells_fullcase[I00{idx}.get()]);
      _idxs_local.push_back(idx);
    }
  }

  void tearDown() {
    IndexingService<3>::getInstance().finalize();
    for (auto cell : _couplingCells_fullcase) {
      if (cell != nullptr)
        delete cell;
    }
    _couplingCells_fullcase.clear();
  }

  void testInsertAccessSize() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container;
    int i = 0;
    for (auto idx : I01()) {
      (void)idx; // Avoid unused variable error
      container << _couplingCells_fullcase[i++];
    }
    CPPUNIT_ASSERT_EQUAL((unsigned int)_couplingCells_fullcase.size(), container.size());
    auto extractedCell = (*container.begin()).first;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(_couplingCells_fullcase[0]->getMacroscopicMass(), extractedCell->getMacroscopicMass(), 1e-16);
    auto extractedIdx = (*container.begin()).second;
    CPPUNIT_ASSERT_EQUAL(extractedIdx, _idxs_fullcase[0]);
  }

  void testIteration() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container;
    int i = 0;
    for (auto idx : I01()) {
      (void)idx; // Avoid unused variable error
      container << _couplingCells_fullcase[i++];
    }
    I01 idx;
    CouplingCell<3>* couplingCell;
    i = 0;
    for (auto pair : container) {
      std::tie(couplingCell, idx) = pair;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_fullcase[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_fullcase[i]);
      i++;
    }
  }

  void testIterationLocal() {
    using namespace coupling::datastructures;
    CellContainer<I03, 3> container(_couplingCells_local);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto pair : container) {
      std::tie(couplingCell, idx) = pair;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_local[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT(idx == _idxs_local[i]);
      i++;
    }
  }

  void testEquality() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container1(_couplingCells_fullcase);
    auto iter1 = container1.begin();
    auto iter2 = container1.begin();
    for (unsigned int i = 0; i < _couplingCells_fullcase.size(); ++i) {
      CPPUNIT_ASSERT(iter1 == iter2);
      iter1++;
      iter2++;
    }
  }

  void testEqualityContent() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container1(_couplingCells_fullcase);
    CellContainer<I01, 3> container2(_couplingCells_fullcase);
    auto iter1 = container1.begin();
    auto iter2 = container2.begin();
    for (unsigned int i = 0; i < _couplingCells_fullcase.size(); ++i) {
      CPPUNIT_ASSERT(*iter1 == *iter2);
      iter1++;
      iter2++;
    }
  }

  void testInequality() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container1(_couplingCells_fullcase);
    CellContainer<I01, 3> container2(_couplingCells_fullcase);
    auto iter1 = container1.begin();
    auto iter2 = container2.begin();
    for (unsigned int i = 0; i < _couplingCells_fullcase.size(); ++i) {
      iter1++;
      CPPUNIT_ASSERT(iter1 != iter2);
      iter2++;
    }
  }

  void testPreIncrement() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container(_couplingCells_fullcase);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto iter = container.begin(); iter != container.end();) {
      std::tie(couplingCell, idx) = *iter;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_fullcase[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_fullcase[i]);
      i++;
      ++iter;
    }
  }

  void testPostIncrement() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container(_couplingCells_fullcase);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto iter = container.begin(); iter != container.end();) {
      std::tie(couplingCell, idx) = *iter++;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_fullcase[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_fullcase[i]);
      i++;
    }
  }

  void testIndexAccess() {
    using namespace coupling::datastructures;
    CellContainer<I01, 3> container(_couplingCells_fullcase);
    I01 idx;
    unsigned int i = 0;
    for (; i < _couplingCells_fullcase.size(); i++) {
      auto idx = _idxs_fullcase[i];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(container[idx]->getMacroscopicMass(), _couplingCells_fullcase[i]->getMacroscopicMass(), 1e-16);
    }
  }

private:
  int _size, _rank;
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells_fullcase;
  std::vector<I01> _idxs_fullcase;
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells_local;
  std::vector<I03> _idxs_local;
};

CPPUNIT_TEST_SUITE_REGISTRATION(CellContainerTest);
