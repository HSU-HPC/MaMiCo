#include "coupling/datastructures/LinkedCellContainer.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSolverInterface.h"
#include "simplemd/LinkedCell.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

using namespace coupling::indexing;

class LinkedCellContainerTest : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(LinkedCellContainerTest);
  CPPUNIT_TEST(testApplyToLocalNonGhostCouplingCellsWithLinkedCells);
  CPPUNIT_TEST(testApplyToLocalGhostCouplingCellsWithLinkedCells);
  CPPUNIT_TEST(testApplyToAllLocalCouplingCellsWithLinkedCells);
  CPPUNIT_TEST(testApplyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells);
  CPPUNIT_TEST(testApplyXLayersOfGlobalNonGhostCellsWithLinkedCells);
  CPPUNIT_TEST_SUITE_END();

private:
  /** test interface; only implements the getLinkedCell method (required in initialisation of macroscopic cells) and
   *  strictly returns its only linked cell ;-) The MD domain is considered to be the unit square/box.
   */
  class TestMDSolverInterface : public coupling::interface::MDSolverInterface<simplemd::LinkedCell, 3> {
  public:
    TestMDSolverInterface() : coupling::interface::MDSolverInterface<simplemd::LinkedCell, 3>() {}
    virtual ~TestMDSolverInterface() {}

    virtual simplemd::LinkedCell& getLinkedCell(const I11& couplingCellIndex, const tarch::la::Vector<3, unsigned int>& linkedCellInCouplingCell,
                                                const tarch::la::Vector<3, unsigned int>& linkedCellsPerCouplingCell) {
      return _linkedcell;
    }

    /** returns the global size of the box-shaped MD domain */
    virtual tarch::la::Vector<3, double> getGlobalMDDomainSize() const { return tarch::la::Vector<3, double>(30.0); }

    /** returns the offset (i.e. lower,left corner) of MD domain */
    virtual tarch::la::Vector<3, double> getGlobalMDDomainOffset() const { return tarch::la::Vector<3, double>(0.0); }

    virtual double getMoleculeMass() const { return 1.0; }
    virtual double getKB() const { return 1.0; }
    virtual double getMoleculeSigma() const { return 1.0; }
    virtual double getMoleculeEpsilon() const { return 1.0; }
    virtual void getInitialVelocity(const tarch::la::Vector<3, double>& meanVelocity, const double& kB, const double& temperature,
                                    tarch::la::Vector<3, double>& initialVelocity) const {}
    virtual void deleteMoleculeFromMDSimulation(const coupling::interface::Molecule<3>& molecule, simplemd::LinkedCell& cell) {}
    virtual void addMoleculeToMDSimulation(const coupling::interface::Molecule<3>& molecule) {}
    virtual void setupPotentialEnergyLandscape(const tarch::la::Vector<3, unsigned int>& indexOfFirstCouplingCell,
                                               const tarch::la::Vector<3, unsigned int>& rangeCouplingCells,
                                               const tarch::la::Vector<3, unsigned int>& linkedCellsPerCouplingCell) {}
    virtual tarch::la::Vector<3, unsigned int> getLinkedCellIndexForMoleculePosition(const tarch::la::Vector<3, double>& position) {
      return tarch::la::Vector<3, unsigned int>(0);
    }
    virtual void calculateForceAndEnergy(coupling::interface::Molecule<3>& molecule) {}
    virtual void synchronizeMoleculesAfterMassModification() {}
    virtual void synchronizeMoleculesAfterMomentumModification() {}
    virtual double getDt() { return 1.0; }
    virtual coupling::interface::MoleculeIterator<simplemd::LinkedCell, 3>* getMoleculeIterator(simplemd::LinkedCell& cell) { return NULL; }

  private:
    simplemd::LinkedCell _linkedcell;
  };

  class LayerChecker {
  public:
    LayerChecker(bool local, std::function<bool(int)> check) : _local{local}, _check{check} {}
    ~LayerChecker() {}
    void beginCellIteration() { _applicationCount = 0; }
    void endCellIteration() {}
    void apply(coupling::datastructures::CouplingCellWithLinkedCells<simplemd::LinkedCell, 3>& cell, I02 index) {
      _applicationCount++;
      CPPUNIT_ASSERT(_check(getLayer(index)));
    }
    int getApplicationCount() { return _applicationCount; }

  private:
    /**
     * Recursive function to find out how far from the domain boundary a cell is
     * Returns 0 for outermost layer of cells in all directions
     * Returns 1 for next-to-outermost layer and so on.
     */
    int getLayer(tarch::la::Vector<3, int> index, tarch::la::Vector<3, unsigned int> size) {
      for (int d = 0; d < 3; d++) {
        CPPUNIT_ASSERT(index[d] >= 0);
        CPPUNIT_ASSERT(index[d] <= static_cast<int>(size[d]) - 1);
        if (index[d] == 0)
          return 0;
        if (index[d] == static_cast<int>(size[d]) - 1)
          return 0;
        CPPUNIT_ASSERT(size[d] > 2);
      }
      size -= tarch::la::Vector<3, unsigned int>{2, 2, 2};
      index -= tarch::la::Vector<3, int>{1, 1, 1};
      return 1 + getLayer(index, size);
    }

    /* Function to find out how far from the domain boundary a cell is
     */
    int getLayer(I02 index) {
      if (_local)
        return getLayer(CellIndex<3, IndexTrait::vector, IndexTrait::local>{index}.get(),
                        CellIndex<3, IndexTrait::vector, IndexTrait::local>::numberCellsInDomain);
      else
        return getLayer(CellIndex<3, IndexTrait::vector>{index}.get(), CellIndex<3, IndexTrait::vector>::numberCellsInDomain);
    }
    int _applicationCount;
    bool _local;
    std::function<bool(int)> _check;
  };

public:
  void setUp() {
    _rank = 0;
    _size = 1;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
#endif

    const tarch::la::Vector<3, unsigned int> globalNumberCells(12);
    tarch::la::Vector<3, unsigned int> numberProcesses{(unsigned int)_size, 1, 1};
    _testInterface = new TestMDSolverInterface();
    IndexingService<3>::getInstance().initWithCells(globalNumberCells, numberProcesses, coupling::paralleltopology::XYZ, 3, (unsigned int)_rank);
    const tarch::la::Vector<3, unsigned int> numberLinkedCellsPerCouplingCell(2);
    _cells = std::make_unique<coupling::datastructures::LinkedCellContainer<simplemd::LinkedCell, 3>>(numberLinkedCellsPerCouplingCell, _testInterface);
  }

  void tearDown() {
    delete _testInterface;
    IndexingService<3>::getInstance().finalize();
  }

  void testApplyToLocalNonGhostCouplingCellsWithLinkedCells() {
    LayerChecker check{true, [](int layer) { return layer > 0; }};
    _cells->applyToLocalNonGhostCouplingCellsWithLinkedCells(check);
    CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 12 * 12 * 12 / _size);
  }

  void testApplyToLocalGhostCouplingCellsWithLinkedCells() {
    LayerChecker check{true, [](int layer) { return layer == 0; }};
    _cells->applyToLocalGhostCouplingCellsWithLinkedCells(check);
    auto numlocal = CellIndex<3, IndexTrait::local>::numberCellsInDomain;
    int numghost = numlocal[0] * numlocal[1] * numlocal[2] - (numlocal[0] - 2) * (numlocal[1] - 2) * (numlocal[2] - 2);
    CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), numghost);
  }

  void testApplyToAllLocalCouplingCellsWithLinkedCells() {
    LayerChecker check{true, [](int layer) { return true; }};
    _cells->applyToAllLocalCouplingCellsWithLinkedCells(check);
    int numlocal = CellIndex<3, IndexTrait::local>::linearNumberCellsInDomain;
    CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), numlocal);
  }

  void testApplyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells() {
    LayerChecker check{false, [](int layer) { return layer == 1; }};
    _cells->applyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells(check);
    auto numcells = CellIndex<3, IndexTrait::noGhost>::numberCellsInDomain;
    int num1stLayer = numcells[0] * numcells[1] * numcells[2] - (numcells[0] - 2) * (numcells[1] - 2) * (numcells[2] - 2);
    if (_size == 1)
      CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), num1stLayer);
    if (_size == 2)
      CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), num1stLayer / 2);
    if (_size == 4) {
      if (_rank == 0 || _rank == 3)
        CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 232);
      if (_rank == 1 || _rank == 2)
        CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 132);
    }
  }

  void testApplyXLayersOfGlobalNonGhostCellsWithLinkedCells() {
    for (int x = 1; x <= 6; x++) {
      LayerChecker check{false, [x](int layer) { return layer > 0 && layer <= x; }};
      _cells->applyXLayersOfGlobalNonGhostCellsWithLinkedCells(check, x);
      auto numcells = CellIndex<3, IndexTrait::noGhost>::numberCellsInDomain;
      int numXLayer = numcells[0] * numcells[1] * numcells[2] - (numcells[0] - 2 * x) * (numcells[1] - 2 * x) * (numcells[2] - 2 * x);
      if (_size == 1)
        CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), numXLayer);
      if (_size == 2)
        CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), numXLayer / 2);
      if (_size == 4) {
        if (x == 1) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 232);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 132);
        }
        if (x == 2) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 368);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 240);
        }
        if (x == 3) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 432);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 324);
        }
        if (x == 4) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 432);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 400);
        }
        if (x == 5) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 432);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 428);
        }
        if (x == 6) {
          if (_rank == 0 || _rank == 3)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 432);
          if (_rank == 1 || _rank == 2)
            CPPUNIT_ASSERT_EQUAL(check.getApplicationCount(), 432);
        }
      }
    }
  }

private:
  int _size, _rank;
  std::unique_ptr<coupling::datastructures::LinkedCellContainer<simplemd::LinkedCell, 3>> _cells;
  TestMDSolverInterface* _testInterface;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LinkedCellContainerTest);