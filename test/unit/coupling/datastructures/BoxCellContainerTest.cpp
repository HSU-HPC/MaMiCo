#include "coupling/datastructures/BoxCellContainer.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <tuple>
#include <vector>

using namespace coupling::indexing;

class BoxCellContainerTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(BoxCellContainerTest);
  CPPUNIT_TEST(testIteration);
  CPPUNIT_TEST(testIteration2);
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

  void testIteration() {
    using namespace coupling::datastructures;
    BoxCellContainer container(*I01::begin(), I01::numberCellsInDomain);
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

    // TODO write more tests
    // TODO modify setUp and tearDown to use other position/size (not indexing domains)

  }

  void testIteration2(){
    I01 lowerBound{1,2,3};
    tarch::la::Vector<3, int> shape{1,2,0};
    BoxCellContainer box{lowerBound, shape};
    for(auto cell : box)
      CPPUNIT_FAIL("Iterated over a cell, but container is empty!");
    tarch::la::Vector<3, int> shape2{1,2,3};
    BoxCellContainer box2{lowerBound, shape2};
    int cnt = 0;
    for(auto cell : box2)
      cnt++
    CPPUNIT_ASSERT_EQUAL(cnt, 6);
  }

private:
  int _size, _rank;
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells_fullcase;
  std::vector<I01> _idxs_fullcase;
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells_local;
  std::vector<I03> _idxs_local;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BoxCellContainerTest);
