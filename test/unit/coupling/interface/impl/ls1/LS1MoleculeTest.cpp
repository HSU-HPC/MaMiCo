#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"
// add simulation in includes! dont rely on include chain!
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/impl/ls1/LS1Molecule.h"

#include "ls1/src/Simulation.h"

#include <sstream>

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1MoleculeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LS1MoleculeTest);
  CPPUNIT_TEST(testPosition);
  CPPUNIT_TEST(testVelocity);
  CPPUNIT_TEST(testForce);
  CPPUNIT_TEST(testGetPotentialEnergy);
  CPPUNIT_TEST(testSetPotentialEnergy);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
#ifdef MARDYN_AUTOPAS
    _ls1ConfigFileName = "../test/unit/coupling/interface/impl/ls1/autopasgridconfig.xml";
#else
    _ls1ConfigFileName = "../test/unit/coupling/interface/impl/ls1/ls1gridconfig.xml";
#endif
    Log::global_log = std::make_unique<Log::Logger>(Log::None);
    Log::global_log->set_mpi_output_root(0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, 0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, 0);
    coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, 0);
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    coupling::interface::LS1StaticCommData::getInstance().setLocalCommunicator(MPI_COMM_WORLD);
#endif
    _testSimulation = new Simulation();
    global_simulation = _testSimulation;
    _testSimulation->disableFinalCheckpoint();
    _testSimulation->readConfigFile(_ls1ConfigFileName);
    //_testSimulation->getDomain()->thermostatOff();
    //_testSimulation->getDomain()->setExplosionHeuristics(false);
    // after this point the mamico plugin exists and is accessible
    _testSimulation->prepare_start();
    _testSimulation->preSimLoopSteps();
  }
  void tearDown() {}
  void testPosition() {
    ::Molecule ls1Molecule;
    coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
    for (int i = 4; i < 8; i++) {
      coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(0, i + 0.5);
      for (int j = 4; j < 8; j++) {
        coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(1, j + 0.5);
        for (int k = 4; k < 8; k++) {
          coupling::interface::LS1StaticCommData::getInstance().setBoxOffsetAtDim(2, k + 0.5);
          tarch::la::Vector<3, double> globalOffset{i + 0.5, j + 0.5, k + 0.5};
          for (int ii = 4; ii < 8; ii++) {
            for (int jj = 4; jj < 8; jj++) {
              for (int kk = 4; kk < 8; kk++) {
                tarch::la::Vector<3, double> position{i + 0.5, j + 0.5, k + 0.5};

                mamicoMolecule.setPosition(position);
                tarch::la::Vector<3, double> storedPosition = mamicoMolecule.getPosition();

                tarch::la::Vector<3, double> internalPosition{ls1Molecule.r(0), ls1Molecule.r(1), ls1Molecule.r(2)};

                std::stringstream storedPosStr;
                storedPosStr << "stored position assertion-> storedPosition: " << storedPosition << " position: " << position << std::endl;
                std::stringstream offsetCalcStr;
                offsetCalcStr << "offset calculation-> storedPosition: " << storedPosition << " internal position: " << internalPosition
                              << " global offset: " << globalOffset << std::endl;

                CPPUNIT_ASSERT_MESSAGE(storedPosStr.str(), storedPosition == position);
                CPPUNIT_ASSERT_MESSAGE(offsetCalcStr.str(), storedPosition == internalPosition + globalOffset);
              }
            }
          }
        }
      }
    }
  }
  void testVelocity() {
    ::Molecule ls1Molecule;
    coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
    for (int i = -4000; i < 4001; i += 100) {
      for (int j = -4000; j < 4001; j += 100) {
        for (int k = -4000; k < 4001; k += 100) {
          tarch::la::Vector<3, double> velocity{i + 0.5, j + 0.5, k + 0.5};

          mamicoMolecule.setVelocity(velocity);
          tarch::la::Vector<3, double> storedVelocity = mamicoMolecule.getVelocity();
          CPPUNIT_ASSERT_MESSAGE("velocity assertion", storedVelocity == velocity);
        }
      }
    }
  }

  void testForce() {
    ::Molecule ls1Molecule;
    coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
    for (int i = -4000; i < 4001; i += 100) {
      for (int j = -4000; j < 4001; j += 100) {
        for (int k = -4000; k < 4001; k += 100) {
          tarch::la::Vector<3, double> force{i + 0.5, j + 0.5, k + 0.5};

          mamicoMolecule.setForce(force);
          tarch::la::Vector<3, double> storedForce = mamicoMolecule.getForce();
          CPPUNIT_ASSERT_MESSAGE("force assertion", storedForce == force);
        }
      }
    }
  }

  void testGetPotentialEnergy() {
    ::Molecule ls1Molecule;
    coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
    ls1::LS1RegionWrapper region;
    for(int i = 0; i < 40; i++) {
      for(int j = 0; j < 40; j++) {
        for(int k = 0; k < 40; k++) {
          tarch::la::Vector<3, double> position{i + 0.5, j + 0.5, k + 0.5};
          mamicoMolecule.setPosition(position);
          const double potFromMolecule = mamicoMolecule.getPotentialEnergy();
          const auto [force, potFromIterator] = region.calculateForceAndPotentialAtPoint(position, true, true, false);
          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Potential assersion", potFromIterator, potFromMolecule, 1e-6);
        }
      }
    }
  }

  void testSetPotentialEnergy() {
    ::Molecule ls1Molecule;
    coupling::interface::LS1Molecule mamicoMolecule(&ls1Molecule);
    CPPUNIT_ASSERT_THROW(mamicoMolecule.setPotentialEnergy(0.0), std::runtime_error);
  }

private:
  Simulation* _testSimulation;
  std::string _ls1ConfigFileName;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1MoleculeTest);