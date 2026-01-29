#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/MoleculeContainer.h"
#include "simplemd/LinkedCell.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/la/Vector.h"

#include <ctime>
#include <array>

template <bool par> class MolPosIncrMapping {
public:
  MolPosIncrMapping() {}
  KOKKOS_FUNCTION void beginMoleculeIteration() {};
  KOKKOS_FUNCTION void endMoleculeIteration() {};
  KOKKOS_FUNCTION void handleMolecule(simplemd::Molecule& molecule) const {
    tarch::la::Vector<MD_DIM, double>& pos = molecule.getPosition();
    for (size_t i = 0; i < MD_DIM; i++) {
      pos[i] += i + 1;
    }
  }
  static const bool IsParallel = par;
};

template <bool par> class CellMolPosIncrMapping {
public:
  CellMolPosIncrMapping() {}
  KOKKOS_FUNCTION void beginCellIteration() {};
  KOKKOS_FUNCTION void endCellIteration() {};
  KOKKOS_FUNCTION void handleCell(simplemd::LinkedCell& cell) const {
    for (auto it = cell.begin(); it != cell.end(); it++) {
      tarch::la::Vector<MD_DIM, double>& pos = it->getPosition();
      for (size_t i = 0; i < MD_DIM; i++) {
        pos[i] += i + 1;
      }
    }
  }
  static const bool IsParallel = par;
};

template <bool par> class CellPairPosMaxMapping {
public:
  CellPairPosMaxMapping() {}
  KOKKOS_FUNCTION void beginCellIteration() {};
  KOKKOS_FUNCTION void endCellIteration() {};
  KOKKOS_FUNCTION void handleCell(simplemd::LinkedCell& cell) const {};
  KOKKOS_FUNCTION void handleCellPair(simplemd::LinkedCell& cell1, simplemd::LinkedCell& cell2, const unsigned int& cellIndex1,
                                      const unsigned int& cellIndex2) const {
    double temp;
    for (auto it1 = cell1.begin(); it1 != cell1.end(); it1++) {
      tarch::la::Vector<MD_DIM, double>& pos1 = it1->getPosition();
      for (auto it2 = cell2.begin(); it2 != cell2.end(); it2++) {
        tarch::la::Vector<MD_DIM, double>& pos2 = it2->getPosition();
        for (size_t i = 0; i < MD_DIM; i++) {
          temp = pos1[i] > pos2[i] ? pos1[i] : pos2[i];
          pos1[i] = temp;
          pos2[i] = temp;
        }
      }
    }
  }
  static const bool IsParallel = par;
};

class MoleculeContainerTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(MoleculeContainerTest);
  CPPUNIT_TEST(testInsertRemove);
  CPPUNIT_TEST(testSort);
  CPPUNIT_TEST(testClearCell);
  CPPUNIT_TEST(testGetMoleculeAt);
  CPPUNIT_TEST(testNumGhostCells);
  CPPUNIT_TEST(testIterationMoleculesSerial);
  CPPUNIT_TEST(testIterationMoleculesParallel);
  CPPUNIT_TEST(testIterationMoleculesAlternate);
  CPPUNIT_TEST(testIterationLinkedCellsSerial);
  CPPUNIT_TEST(testIterationLinkedCellsParallel);
  CPPUNIT_TEST(testIterationLinkedCellPairsSerial);
  CPPUNIT_TEST(testIterationLinkedCellPairsParallel);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    tarch::la::Vector<MD_DIM, double> numCells(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      numCells[i] = _numCellsIf3D[i];
    }
    const tarch::la::Vector<MD_DIM, double> domainOffset(0);
    const tarch::la::Vector<MD_DIM, double> meshWidth(1);
    const tarch::la::Vector<MD_DIM, unsigned int> numberProcesses(1);
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> boundary(simplemd::BoundaryType::PERIODIC_BOUNDARY);
    simplemd::services::ParallelTopologyService parallelTopologyService(numCells, domainOffset, meshWidth, numberProcesses, boundary);
    _moleculeContainer = new simplemd::MoleculeContainer(parallelTopologyService, 20);
  }

  void tearDown() {
    if (_moleculeContainer != nullptr) {
      delete _moleculeContainer;
      _moleculeContainer = nullptr;
    }
  }

  void testInsertRemove() {
    tarch::la::Vector<MD_DIM, double> position, velocity(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      position[i] = 500;
    }
    simplemd::Molecule trialMolecule(position, velocity);
    _moleculeContainer->insert(0, trialMolecule);
    auto extractMolecule = _moleculeContainer->getMoleculeAt(0, 0);
    for (size_t i = 0; i < MD_DIM; i++) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule.getPosition()[i], extractMolecule.getPosition()[i], 1e-6);
    }
    _moleculeContainer->remove(0, 0);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 0);

    for (size_t i = 0; i < MD_DIM; i++) {
      position[i] = 0.5;
    }
    simplemd::Molecule trialMolecule2(position, velocity);
    _moleculeContainer->insert(trialMolecule2);
    simplemd::Molecule extractMolecule2;
    int index;
#if (MD_DIM == 1)
    // molecule idx is [1] thus cell 1
    index = 1;
#endif
#if (MD_DIM == 2)
    // molecule idx is [1,1]
    index = _numCellsIf3D[0] + 2 + 1;
#endif
#if (MD_DIM == 3)
    // molecule idx is [1,1,1]
    index = (_numCellsIf3D[1] + 2) * (_numCellsIf3D[0] + 2) + _numCellsIf3D[0] + 2 + 1;
#endif

    extractMolecule2 = _moleculeContainer->getMoleculeAt(index, 0);
    for (size_t i = 0; i < MD_DIM; i++) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule2.getPosition()[i], extractMolecule2.getPosition()[i], 1e-6);
    }
  }

  void testSort() {
    const unsigned int numCellsForTest = 10;
    tarch::la::Vector<MD_DIM, double> position, velocity(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      position[i] = 0.5;
    }
    std::array<simplemd::Molecule*, numCellsForTest> molecules;
    for (size_t i = 0; i < numCellsForTest; i++) {
      molecules[i] = new simplemd::Molecule(position, velocity);
      (*_moleculeContainer)[i].insert(*molecules[i]);
    }
    _moleculeContainer->sort();
    // all molecules should be in the same cell since all have same position
    int index;
#if (MD_DIM == 1)
    // molecule idx is [1] thus cell 1
    index = 1;
#endif
#if (MD_DIM == 2)
    // molecule idx is [1,1]
    index = _numCellsIf3D[0] + 2 + 1;
#endif
#if (MD_DIM == 3)
    // molecule idx is [1,1,1]
    index = (_numCellsIf3D[1] + 2) * (_numCellsIf3D[0] + 2) + _numCellsIf3D[0] + 2 + 1;
#endif
    CPPUNIT_ASSERT_EQUAL(numCellsForTest, (*_moleculeContainer)[index].numMolecules());
    // cleanup
    for (size_t i = 0; i < numCellsForTest; i++) {
      _moleculeContainer->clearLinkedCell(i);
      if (molecules[i] != nullptr) {
        delete molecules[i];
        molecules[i] = nullptr;
      }
    }
  }

  void testClearCell() {
    tarch::la::Vector<MD_DIM, double> position1, position2, position3, velocity(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      position1[i] = 500;
      position2[i] = 600;
      position3[i] = 700;
    }
    simplemd::Molecule mol1(position1, velocity), mol2(position2, velocity), mol3(position3, velocity);
    _moleculeContainer->insert(1, mol1);
    _moleculeContainer->insert(1, mol2);
    _moleculeContainer->insert(2, mol3);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 3);
    _moleculeContainer->clearLinkedCell(2);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 2);
    _moleculeContainer->clearLinkedCell(1);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 0);
  }

  void testGetMoleculeAt() {
    tarch::la::Vector<MD_DIM, double> position1, position2, position3, velocity(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      position1[i] = 500;
      position2[i] = 600;
      position3[i] = 700;
    }
    simplemd::Molecule mol1(position1, velocity), mol2(position2, velocity), mol3(position3, velocity);
    _moleculeContainer->insert(1, mol1);
    _moleculeContainer->insert(1, mol2);
    _moleculeContainer->insert(2, mol3);
    for (size_t i = 0; i < MD_DIM; i++) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(mol1.getPosition()[i], _moleculeContainer->getMoleculeAt(1, 0).getPosition()[i], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(mol2.getPosition()[i], _moleculeContainer->getMoleculeAt(1, 1).getPosition()[i], 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(mol3.getPosition()[i], _moleculeContainer->getMoleculeAt(2, 0).getPosition()[i], 1e-6);
    }
    // cleanup
    _moleculeContainer->clearLinkedCell(1);
    _moleculeContainer->clearLinkedCell(2);
  }

  void testNumGhostCells() {
    srand(time(NULL));

    const int numTests = 20;
    tarch::la::Vector<MD_DIM, double> numCellsTestInput;

    const tarch::la::Vector<MD_DIM, double> domainOffset(0);
    const tarch::la::Vector<MD_DIM, double> meshWidth(1);
    const tarch::la::Vector<MD_DIM, unsigned int> numberProcesses(1);
    unsigned int numInnerCells, numTotalCells;

    for (int i = 0; i < numTests; i++) {
      numInnerCells = 1;
      numTotalCells = 1;
      for (int j = 0; j < MD_DIM; j++) {
        numCellsTestInput[j] = static_cast<double>(rand() % 50 + 1);
        numInnerCells *= static_cast<int>(numCellsTestInput[j]);
        numTotalCells *= static_cast<int>(numCellsTestInput[j] + 2);
      }
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> boundary(simplemd::BoundaryType::PERIODIC_BOUNDARY);
      simplemd::services::ParallelTopologyService parallelTopologyService(numCellsTestInput, domainOffset, meshWidth, numberProcesses, boundary);
      simplemd::MoleculeContainer moleculeContainer(parallelTopologyService, 20);
      CPPUNIT_ASSERT(numTotalCells == moleculeContainer.getLocalNumberOfCellsScalarWithGhost());
      CPPUNIT_ASSERT(numInnerCells == vectorToScalar(moleculeContainer.getLocalNumberOfCells()));
    }
  }

  void testIterationMoleculesSerial() { consolidatedMolIter<false>(); }
  void testIterationMoleculesParallel() { consolidatedMolIter<true>(); }

  template <bool par> void consolidatedMolIter() {
    MolPosIncrMapping<par> mapping;
    const unsigned int numCellsForTest = 10;
    const unsigned int numMolsForTestPerCell = 5;
    // generate particle positions
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<std::array<tarch::la::Vector<MD_DIM, double>, numMolsForTestPerCell>, numCellsForTest> positions;
    std::array<std::array<simplemd::Molecule*, numMolsForTestPerCell>, numCellsForTest> molecules;
    int ctr = 5; // arbitrary value
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          positions[i][j][k] = ctr++;
        }
        molecules[i][j] = new simplemd::Molecule(positions[i][j], velocity);
        (*_moleculeContainer)[i].insert(*molecules[i][j]);
      }
    }

    // use mapping
    _moleculeContainer->iterateMolecules(mapping);

    // check if successful
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(positions[i][j][k] + k + 1, _moleculeContainer->getMoleculeAt(i, j).getConstPosition()[0], 1e6);
        }
      }
    }

    // cleanup
    for (size_t i = 0; i < numCellsForTest; i++) {
      _moleculeContainer->clearLinkedCell(i);
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        if (molecules[i][j] != nullptr) {
          delete molecules[i][j];
          molecules[i][j] = nullptr;
        }
      }
    }
  }

  void testIterationMoleculesAlternate() {
    MolPosIncrMapping<true> mappingPar;
    MolPosIncrMapping<false> mappingSer;
    const unsigned int numCellsForTest = 10;
    const unsigned int numMolsForTestPerCell = 5;
    // generate particle positions
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<std::array<tarch::la::Vector<MD_DIM, double>, numMolsForTestPerCell>, numCellsForTest> positions;
    std::array<std::array<simplemd::Molecule*, numMolsForTestPerCell>, numCellsForTest> molecules;
    int ctr = 5; // arbitrary value
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          positions[i][j][k] = ctr++;
        }
        molecules[i][j] = new simplemd::Molecule(positions[i][j], velocity);
        (*_moleculeContainer)[i].insert(*molecules[i][j]);
      }
    }

    // use mapping
    _moleculeContainer->iterateMolecules(mappingSer);
    _moleculeContainer->iterateMolecules(mappingPar);
    _moleculeContainer->iterateMolecules(mappingSer);
    _moleculeContainer->iterateMolecules(mappingPar);

    // check if successful
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(positions[i][j][k] + 4 * (k + 1), _moleculeContainer->getMoleculeAt(i, j).getConstPosition()[0], 1e6);
        }
      }
    }

    // cleanup
    for (size_t i = 0; i < numCellsForTest; i++) {
      _moleculeContainer->clearLinkedCell(i);
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        if (molecules[i][j] != nullptr) {
          delete molecules[i][j];
          molecules[i][j] = nullptr;
        }
      }
    }
  }

  void testIterationLinkedCellsSerial() { consolidatedCellIter<false>(); }
  void testIterationLinkedCellsParallel() { consolidatedCellIter<true>(); }

  template <bool par> void consolidatedCellIter() {
    CellMolPosIncrMapping<par> mapping;
    const unsigned int numCellsForTest = 10;
    const unsigned int numMolsForTestPerCell = 5;
    // generate particle positions
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<std::array<tarch::la::Vector<MD_DIM, double>, numMolsForTestPerCell>, numCellsForTest> positions;
    std::array<std::array<simplemd::Molecule*, numMolsForTestPerCell>, numCellsForTest> molecules;
    int ctr = 5; // arbitrary value
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          positions[i][j][k] = ctr++;
        }
        molecules[i][j] = new simplemd::Molecule(positions[i][j], velocity);
        (*_moleculeContainer)[i].insert(*molecules[i][j]);
      }
    }

    // use mapping
    _moleculeContainer->iterateCells(mapping);

    // check if successful
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(positions[i][j][k] + k + 1, _moleculeContainer->getMoleculeAt(i, j).getConstPosition()[0], 1e6);
        }
      }
    }

    // cleanup
    for (size_t i = 0; i < numCellsForTest; i++) {
      _moleculeContainer->clearLinkedCell(i);
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        if (molecules[i][j] != nullptr) {
          delete molecules[i][j];
          molecules[i][j] = nullptr;
        }
      }
    }
  }

  void testIterationLinkedCellPairsSerial() { consolidatedCellPairIter<false>(); };
  void testIterationLinkedCellPairsParallel() { consolidatedCellPairIter<true>(); };

  template <bool par> void consolidatedCellPairIter() {
    CellPairPosMaxMapping<par> mapping;
    // cells 0 and 1 neighbours in all MD_SIM values
    const unsigned int numCellsForTest = 2;
    const unsigned int numMolsForTestPerCell = 15;
    // generate particle positions
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<std::array<tarch::la::Vector<MD_DIM, double>, numMolsForTestPerCell>, numCellsForTest> positions;
    std::array<std::array<simplemd::Molecule*, numMolsForTestPerCell>, numCellsForTest> molecules;
    tarch::la::Vector<MD_DIM, double> maxPos(0);
    int ctr = 5; // arbitrary value
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          positions[i][j][k] = ctr++;
          maxPos[k] = positions[i][j][k] > maxPos[k] ? positions[i][j][k] : maxPos[k]; // always true, but kept for clarity
        }
        molecules[i][j] = new simplemd::Molecule(positions[i][j], velocity);
        (*_moleculeContainer)[i].insert(*molecules[i][j]);
      }
    }

    // use mapping
    _moleculeContainer->iterateCellPairs(mapping);

    // check if successful
    for (size_t i = 0; i < numCellsForTest; i++) {
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        for (size_t k = 0; k < MD_DIM; k++) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(maxPos[k], _moleculeContainer->getMoleculeAt(i, j).getConstPosition()[0], 1e6);
        }
      }
    }

    // cleanup
    for (size_t i = 0; i < numCellsForTest; i++) {
      _moleculeContainer->clearLinkedCell(i);
      for (size_t j = 0; j < numMolsForTestPerCell; j++) {
        if (molecules[i][j] != nullptr) {
          delete molecules[i][j];
          molecules[i][j] = nullptr;
        }
      }
    }
  }

private:
  unsigned int vectorToScalar(tarch::la::Vector<MD_DIM, unsigned int> vector) const {
    unsigned int toRet = 1;
    for (int i = 0; i < MD_DIM; i++) {
      toRet *= vector[i];
    }
    return toRet;
  }

  // use for persistent tests
  const tarch::la::Vector<3, double> _numCellsIf3D = {100, 60, 50};
  simplemd::MoleculeContainer* _moleculeContainer;
};

CPPUNIT_TEST_SUITE_REGISTRATION(MoleculeContainerTest);