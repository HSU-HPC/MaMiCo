#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/MoleculeContainer.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/la/Vector.h"
#include "simplemd/LinkedCell.h"

#include <array>

class LinkedCellTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LinkedCellTest);
  CPPUNIT_TEST(testInsertRemove);
  CPPUNIT_TEST(testClearCell);
  CPPUNIT_TEST(testIteratorIncrement);
  CPPUNIT_TEST(testIteratorDecrement);
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
    (*_moleculeContainer)[0].clear();
    (*_moleculeContainer)[0].insert(trialMolecule);
    auto extractMolecule = _moleculeContainer->getMoleculeAt(0, 0);
    for (size_t i = 0; i < MD_DIM; i++) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule.getPosition()[i], extractMolecule.getPosition()[i], 1e-6);
    }
    (*_moleculeContainer)[0].remove(0);
    CPPUNIT_ASSERT((*_moleculeContainer)[0].numMolecules() == 0);
  }

  void testClearCell() {
    tarch::la::Vector<MD_DIM, double> position1, position2, position3, velocity(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      position1[i] = 500;
      position2[i] = 600;
      position3[i] = 700;
    }
    simplemd::Molecule mol1(position1, velocity), mol2(position2, velocity), mol3(position3, velocity);
    (*_moleculeContainer)[1].insert(mol1);
    (*_moleculeContainer)[1].insert(mol2);
    (*_moleculeContainer)[2].insert(mol3);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 3);
    CPPUNIT_ASSERT((*_moleculeContainer)[1].numMolecules() == 2);
    CPPUNIT_ASSERT((*_moleculeContainer)[2].numMolecules() == 1);
    (*_moleculeContainer)[2].clear();
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 2);
    CPPUNIT_ASSERT((*_moleculeContainer)[2].numMolecules() == 0);
    (*_moleculeContainer)[1].clear();
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 0);
    CPPUNIT_ASSERT((*_moleculeContainer)[1].numMolecules() == 0);
  }

  void testIteratorIncrement() {
    const unsigned int numMols = 5;
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<tarch::la::Vector<MD_DIM, double>, numMols> positions;
    std::array<simplemd::Molecule*, numMols> molecules;
    for (unsigned int i = 0; i < positions.size(); i++) {
      positions[i][0] = i + 1;
      molecules[i] = new simplemd::Molecule(positions[i], velocity);
      (*_moleculeContainer)[0].insert(*molecules[i]);
    }

    // prefix
    double pos = 1;
    for (auto it = (*_moleculeContainer)[0].begin(); it != (*_moleculeContainer)[0].end(); ++it) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
      unsigned int idx = it.getIndex();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);
      pos++;
    }

    // postfix
    pos = 1;
    for (auto it = (*_moleculeContainer)[0].begin(); it != (*_moleculeContainer)[0].end(); it++) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
      unsigned int idx = it.getIndex();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);
      pos++;
    }

    // cleanup
    (*_moleculeContainer)[0].clear();
    for (unsigned int i = 0; i < positions.size(); i++) {
      if (molecules[i] != nullptr) {
        delete molecules[i];
        molecules[i] = nullptr;
      }
    }
  }

  void testIteratorDecrement() {
    const unsigned int numMols = 5;
    tarch::la::Vector<MD_DIM, double> velocity(0);
    std::array<tarch::la::Vector<MD_DIM, double>, numMols> positions;
    std::array<simplemd::Molecule*, numMols> molecules;
    for (unsigned int i = 0; i < positions.size(); i++) {
      positions[i][0] = i + 1;
      molecules[i] = new simplemd::Molecule(positions[i], velocity);
      (*_moleculeContainer)[0].insert(*molecules[i]);
    }
    // prefix
    double pos = numMols;
    auto it = (*_moleculeContainer)[0].end();
    it--;
    for (; it != (*_moleculeContainer)[0].begin(); --it) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
      unsigned int idx = it.getIndex();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);
      pos--;
    }
    // check first item
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
    unsigned int idx = it.getIndex();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);
    // postfix
    it = (*_moleculeContainer)[0].end();
    it--;
    pos = numMols;
    for (; it != (*_moleculeContainer)[0].begin(); it--) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
      unsigned int idx = it.getIndex();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);
      pos--;
    }
    // check first item
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], pos, 1e-6);
    idx = it.getIndex();
    CPPUNIT_ASSERT_DOUBLES_EQUAL(it->getPosition()[0], (_moleculeContainer->getMoleculeAt(0, idx)).getPosition()[0], 1e-6);

    // cleanup
    (*_moleculeContainer)[0].clear();
    for (unsigned int i = 0; i < positions.size(); i++) {
      if (molecules[i] != nullptr) {
        delete molecules[i];
        molecules[i] = nullptr;
      }
    }
  }

private:
  int vectorToScalar(tarch::la::Vector<MD_DIM, unsigned int> vector) const {
    int toRet = 1;
    for (size_t i = 0; i < MD_DIM; i++) {
      toRet *= vector[i];
    }
    return toRet;
  }

  // use for persistent tests
  const tarch::la::Vector<3, double> _numCellsIf3D = {100, 60, 50};
  simplemd::MoleculeContainer* _moleculeContainer;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LinkedCellTest);