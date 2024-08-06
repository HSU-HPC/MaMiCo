#include "coupling/datastructures/FlexibleCellContainer.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <tuple>
#include <vector>

using namespace coupling::indexing;

class FlexibleCellContainerTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(FlexibleCellContainerTest);
  CPPUNIT_TEST(testInsertAndAccess);
  CPPUNIT_TEST(testMassInsert);
  CPPUNIT_TEST(testSize);
  CPPUNIT_TEST(testIteration);
  CPPUNIT_TEST(testEquality);
  CPPUNIT_TEST(testInequality);
  CPPUNIT_TEST(testPreIncrement);
  CPPUNIT_TEST(testPostIncrement);
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
    IndexingService<3>::getInstance().initWithCells(12, numberProcesses, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
    int i = 0;
    for (auto idx : I01()) {
      if (i > 9)
        break;
      _idxs_10case.push_back(idx);
      auto tmp = new CouplingCell<3>;
      tmp->setMacroscopicMass(i++);
      _couplingCells_10case.push_back(tmp);
    }
  }

  void tearDown() {
    IndexingService<3>::getInstance().finalize();
    for (auto cell : _couplingCells_10case) {
      if (cell != nullptr)
        delete cell;
    }
    _couplingCells_10case.clear();
    _idxs_10case.clear();
  }

  void testInsertAndAccess() {
    using namespace coupling::datastructures;
    I01 idx;
    CouplingCell<3>* couplingCell;
    couplingCell->setMacroscopicMass(50.0);
    FlexibleCellContainer<3> container;
    container << std::make_pair(couplingCell, idx);
    CPPUNIT_ASSERT_EQUAL(1, container.size());
    auto extractedCell = (*container.begin()).first;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(50.0, extractedCell->getMacroscopicMass(), 1e-16);
  }

  void testMassInsert() {
    using namespace coupling::datastructures;
    I01 idx;
    CouplingCell<3>* couplingCell;
    FlexibleCellContainer<3> container;
    int i = 0;
    for (; i < _couplingCells_10case.size(); i++) {
      container << std::make_pair(_couplingCells_10case[i], _idxs_10case[i]);
    }
    i = 0;
    for (auto pair : container) {
      std::tie(couplingCell, idx) = pair;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_10case[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_10case[i]);
      i++;
    }
  }

  void testSize() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container(_couplingCells_10case, _idxs_10case);
    int size = _couplingCells_10case.size();
    CPPUNIT_ASSERT_EQUAL(size, container.size());
  }

  void testIteration() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container(_couplingCells_10case, _idxs_10case);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto pair : container) {
      std::tie(couplingCell, idx) = pair;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_10case[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_10case[i]);
      i++;
    }
  }

  void testEquality() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container1(_couplingCells_10case, _idxs_10case);
    FlexibleCellContainer<3> container2(_couplingCells_10case, _idxs_10case);
    auto iter1 = container1.begin();
    auto iter2 = container2.begin();
    for (int i = 0; i < _couplingCells_10case.size(); ++i) {
      CPPUNIT_ASSERT(iter1 == iter2);
      iter1++;
      iter2++;
    }
  }

  void testInequality() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container1(_couplingCells_10case, _idxs_10case);
    FlexibleCellContainer<3> container2(_couplingCells_10case, _idxs_10case);
    auto iter1 = container1.begin();
    auto iter2 = container2.begin();
    for (int i = 0; i < _couplingCells_10case.size(); ++i) {
      iter1++;
      CPPUNIT_ASSERT(iter1 != iter2);
      iter2++;
    }
  }

  void testPreIncrement() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container(_couplingCells_10case, _idxs_10case);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto iter = container.begin(); iter != container.end();) {
      std::tie(couplingCell, idx) = *iter;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_10case[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_10case[i]);
      i++;
      ++iter;
    }
  }

  void testPostIncrement() {
    using namespace coupling::datastructures;
    FlexibleCellContainer<3> container(_couplingCells_10case, _idxs_10case);
    I01 idx;
    CouplingCell<3>* couplingCell;
    int i = 0;
    for (auto iter = container.begin(); iter != container.end();) {
      std::tie(couplingCell, idx) = *iter++;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(couplingCell->getMacroscopicMass(), _couplingCells_10case[i]->getMacroscopicMass(), 1e-16);
      CPPUNIT_ASSERT_EQUAL(idx, _idxs_10case[i]);
      i++;
    }
  }

private:
  int _size, _rank;
  std::vector<I01> _idxs_10case;
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells_10case;
};

CPPUNIT_TEST_SUITE_REGISTRATION(FlexibleCellContainerTest);