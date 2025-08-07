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
  void setUp() {}
  void tearDown() {}
  void testInsertRemove() {}
  void testSort() {}
  void testClearCell() {}
  void testGetMoleculeAt() {}
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
      simplemd::services::ParallelTopologyService parallelTopologyService(numCellsTestInput, domainOffset, meshWidth, numberProcesses,
                                                                          simplemd::BoundaryType::PERIODIC_BOUNDARY);
      simplemd::MoleculeContainer moleculeContainer(parallelTopologyService, 20);
      CPPUNIT_ASSERT((numTotalCells == moleculeContainer.getLocalNumberOfCellsScalarWithGhost()));
      CPPUNIT_ASSERT(numInnerCells == vectorToScalar(moleculeContainer.getLocalNumberOfCells()));
    }
  }

  void testIterationMoleculesSerial() {};
  void testIterationMoleculesParallel() {};
  void testIterationLinkedCellsSerial() {};
  void testIterationLinkedCellsParallel() {};
  void testIterationLinkedCellPairsSerial() {};
  void testIterationLinkedCellPairsParallel() {};

  int vectorToScalar(tarch::la::Vector<MD_DIM, unsigned int> vector) const {
    int toRet = 1;
    for (size_t i = 0; i < MD_DIM; i++) {
      toRet *= vector[i];
    }
    return toRet;
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(MoleculeContainerTest);