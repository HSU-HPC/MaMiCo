#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "simplemd/MoleculeContainer.h"
#include "simplemd/services/ParallelTopologyService.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/la/Vector.h"

#include <sstream>
#include <ctime>
#include <array>

class MoleculeContainerTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(MoleculeContainerTest);
  CPPUNIT_TEST(testInsertRemove);
  CPPUNIT_TEST(testSort);
  CPPUNIT_TEST(testClearCell);
  CPPUNIT_TEST(testGetMoleculeAt);
  CPPUNIT_TEST(testNumGhostCells);
  CPPUNIT_TEST(testIterationMoleculesSerial);
  CPPUNIT_TEST(testIterationMoleculesParallel);
  CPPUNIT_TEST(testIterationLinkedCellsSerial);
  CPPUNIT_TEST(testIterationLinkedCellsParallel);
  CPPUNIT_TEST(testIterationLinkedCellPairsSerial);
  CPPUNIT_TEST(testIterationLinkedCellPairsParallel);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    const tarch::la::Vector<3, double> numCellsIf3D = {100,60,50};
    tarch::la::Vector<MD_DIM, double> numCells(0);
    for (size_t i = 0; i < MD_DIM; i++) {
      numCells[i] = numCellsIf3D[i];
    }
    const tarch::la::Vector<MD_DIM, double> domainOffset(0);
    const tarch::la::Vector<MD_DIM, double> meshWidth(1);
    const tarch::la::Vector<MD_DIM, unsigned int> numberProcesses(1);
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> boundary(simplemd::BoundaryType::PERIODIC_BOUNDARY);
    simplemd::services::ParallelTopologyService parallelTopologyService(numCells, domainOffset, meshWidth, numberProcesses,
                                                                          boundary);
    _moleculeContainer = new simplemd::MoleculeContainer(parallelTopologyService, 20);
  }
  void tearDown() {
    if(_moleculeContainer != nullptr) {
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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule.getPosition()[0], extractMolecule.getPosition()[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule.getPosition()[1], extractMolecule.getPosition()[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(trialMolecule.getPosition()[2], extractMolecule.getPosition()[2], 1e-6);
    _moleculeContainer->remove(0,0);
    CPPUNIT_ASSERT(_moleculeContainer->getLocalNumberOfMoleculesWithGhost() == 0);
  }
  void testSort() {

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
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol1.getPosition()[0], _moleculeContainer->getMoleculeAt(1,0).getPosition()[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol1.getPosition()[1], _moleculeContainer->getMoleculeAt(1,0).getPosition()[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol1.getPosition()[2], _moleculeContainer->getMoleculeAt(1,0).getPosition()[2], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol2.getPosition()[0], _moleculeContainer->getMoleculeAt(1,1).getPosition()[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol2.getPosition()[1], _moleculeContainer->getMoleculeAt(1,1).getPosition()[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol2.getPosition()[2], _moleculeContainer->getMoleculeAt(1,1).getPosition()[2], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol3.getPosition()[0], _moleculeContainer->getMoleculeAt(2,0).getPosition()[0], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol3.getPosition()[1], _moleculeContainer->getMoleculeAt(2,0).getPosition()[1], 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(mol3.getPosition()[2], _moleculeContainer->getMoleculeAt(2,0).getPosition()[2], 1e-6);
    // cleanup
    _moleculeContainer->clearLinkedCell(1);
    _moleculeContainer->clearLinkedCell(2);
  }
  void testNumGhostCells() {
    srand(time(NULL));
    std::stringstream outputMessage;

    const int numTests = 20;
    tarch::la::Vector<MD_DIM, double> numCellsTestInput;

    const tarch::la::Vector<MD_DIM, double> domainOffset(0);
    const tarch::la::Vector<MD_DIM, double> meshWidth(1);
    const tarch::la::Vector<MD_DIM, unsigned int> numberProcesses(1);
    int numInnerCells, numTotalCells;

    for (int i = 0; i < numTests; i++) {
      outputMessage.clear();
      numInnerCells = 1;
      numTotalCells = 1;
      for (int j = 0; j < MD_DIM; j++) {
        numCellsTestInput[j] = static_cast<double>(rand() % 50 + 1);
        numInnerCells *= static_cast<int>(numCellsTestInput[j]);
        numTotalCells *= static_cast<int>(numCellsTestInput[j] + 2);
      }
      const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> boundary(simplemd::BoundaryType::PERIODIC_BOUNDARY);
      simplemd::services::ParallelTopologyService parallelTopologyService(numCellsTestInput, domainOffset, meshWidth, numberProcesses,
                                                                          boundary);
      simplemd::MoleculeContainer moleculeContainer(parallelTopologyService, 20);
      CPPUNIT_ASSERT(numTotalCells == moleculeContainer.getLocalNumberOfCellsScalarWithGhost());
      CPPUNIT_ASSERT(numInnerCells == vectorToScalar(moleculeContainer.getLocalNumberOfCells()));
    }
  }

  void testIterationMoleculesSerial() {};
  void testIterationMoleculesParallel() {};
  void testIterationLinkedCellsSerial() {};
  void testIterationLinkedCellsParallel() {};
  void testIterationLinkedCellPairsSerial() {};
  void testIterationLinkedCellPairsParallel() {};

private:
  int vectorToScalar(tarch::la::Vector<MD_DIM, unsigned int> vector) const {
    int toRet = 1;
    for (size_t i = 0; i < MD_DIM; i++) {
      toRet *= vector[i];
    }
    return toRet;
  }

  //use for persistent tests
  simplemd::MoleculeContainer* _moleculeContainer;
};

CPPUNIT_TEST_SUITE_REGISTRATION(MoleculeContainerTest);