#include "coupling/datastructures/FlexibleCellContainer.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <vector>
#include <tuple>

using namespace coupling::indexing;

class FlexibleCellContainerTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(FlexibleCellContainerTest);
  CPPUNIT_TEST(test);
  CPPUNIT_TEST_SUITE_END();

public:

  void setUp() {
    _rank = 0;
    _size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif
    tarch::la::Vector<3, unsigned int> numberProcesses{1, 1, 1};
    if (_size == 4) numberProcesses = {2, 2, 1};
    IndexingService<3>::getInstance().initWithCells(4, numberProcesses, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
  }

  void tearDown() { IndexingService<3>::getInstance().finalize(); }

  void test() {
    using namespace coupling::datastructures;
    std::vector<I00*> idxs;
    std::vector<CouplingCell<3>*> couplingCells;
    for (auto idx : I00()) {
      // idxs.push_back(&(*idx));
      std::cout << idx << std::endl;
      // couplingCells.push_back(new CouplingCell<3>());
    }
    // FlexibleCellContainer<3> container{couplingCells, idxs};
    // CouplingCell<3>* couplingCell;
    // I00* idx; 
    // for (auto pair : container) {
    //   std::tie(couplingCell, idx) = pair;
    //   std::cout << "idx:" << idx << std::endl;
    // }
    // for (auto idx: idxs) {
    //   if (idx == nullptr) delete idx;
    // }
    // idxs.clear();
    // for (auto couplingCell: couplingCells) {
    //   if (couplingCell == nullptr) delete couplingCell;
    // }
    // couplingCells.clear();
  }

private:
  int _size, _rank;
};

CPPUNIT_TEST_SUITE_REGISTRATION(FlexibleCellContainerTest);