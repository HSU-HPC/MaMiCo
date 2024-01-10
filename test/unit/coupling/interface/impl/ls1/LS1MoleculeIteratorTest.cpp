#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"

#include "coupling/interface/impl/ls1/LS1Molecule.h"
#include "coupling/interface/impl/ls1/LS1MoleculeIterator.h"

#include <sstream>

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

#define MY_LINKEDCELL ls1::LS1RegionWrapper

class LS1MoleculeIteratorTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LS1MoleculeIteratorTest);
  CPPUNIT_TEST(testIteration);
  CPPUNIT_TEST(testReset);
  CPPUNIT_TEST(testGet);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    _ls1ConfigFileName = "../test/unit/coupling/interface/impl/ls1/ls1fourparticleconfig.xml";
    global_log = new Log::Logger(Log::None);
    global_log->set_mpi_output_root(0);
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
  void tearDown() {
    if (_testSimulation != nullptr) {
      delete _testSimulation;
      _testSimulation = nullptr;
    }
  }
  void testIteration() {
    ls1::LS1RegionWrapper fullRegion(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
    coupling::interface::LS1MoleculeIterator moleculeIterator(fullRegion);
    long unsigned int particleCount = 0;
    moleculeIterator.begin();
    while (moleculeIterator.continueIteration()) {
      particleCount++;
      moleculeIterator.next();
    }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    long unsigned int globalCount;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&particleCount, &globalCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    particleCount = globalCount;
#endif
    std::stringstream ss;
    ss << "global particle count: " << particleCount << " particles in simulation: " << _testSimulation->getTotalNumberOfMolecules();
    CPPUNIT_ASSERT_MESSAGE(ss.str(), particleCount == _testSimulation->getTotalNumberOfMolecules());

    // reset
    particleCount = 0;
    moleculeIterator.begin();
    while (moleculeIterator.continueIteration()) {
      particleCount++;
      moleculeIterator.next();
    }

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    globalCount = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&particleCount, &globalCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    particleCount = globalCount;
#endif
    std::stringstream ss2;
    ss2 << "global particle count after reset: " << particleCount << " particles in simulation: " << _testSimulation->getTotalNumberOfMolecules();
    CPPUNIT_ASSERT_MESSAGE(ss2.str(), particleCount == _testSimulation->getTotalNumberOfMolecules());
  }
  void testReset() {
    ls1::LS1RegionWrapper fullRegion(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
    coupling::interface::LS1MoleculeIterator moleculeIterator(fullRegion);
    moleculeIterator.begin();
    if (!moleculeIterator.continueIteration())
      return; // no particles in subdomain

    // grab first particle
    const coupling::interface::Molecule<3>& firstMolecule(moleculeIterator.getConst());

    // iterate through and reset
    while (moleculeIterator.continueIteration()) {
      moleculeIterator.next();
    }

    // reset and check first particle
    moleculeIterator.begin();
    CPPUNIT_ASSERT(firstMolecule.getPosition() == moleculeIterator.getConst().getPosition() &&
                   firstMolecule.getForce() == moleculeIterator.getConst().getForce() &&
                   firstMolecule.getVelocity() == moleculeIterator.getConst().getVelocity());
  }
  void testGet() {
    ls1::LS1RegionWrapper fullRegion(_testSimulation->getEnsemble()->domain()->rmin(), _testSimulation->getEnsemble()->domain()->rmax(), _testSimulation);
    coupling::interface::LS1MoleculeIterator moleculeIterator(fullRegion);
    moleculeIterator.begin();
    if (!moleculeIterator.continueIteration())
      return; // no particles in subdomain

    // grab first particle nonconst
    coupling::interface::Molecule<3>& firstMolecule(moleculeIterator.get());
    firstMolecule.setPosition({0, 0, 0});
    firstMolecule.setForce({0, 0, 0});
    firstMolecule.setVelocity({20, 20, 20});

    // iterate through and reset
    while (moleculeIterator.continueIteration()) {
      moleculeIterator.next();
    }

    // reset and check first particle
    moleculeIterator.begin();
    CPPUNIT_ASSERT(firstMolecule.getPosition() == moleculeIterator.getConst().getPosition() &&
                   firstMolecule.getForce() == moleculeIterator.getConst().getForce() &&
                   firstMolecule.getVelocity() == moleculeIterator.getConst().getVelocity());
  }

private:
  Simulation* _testSimulation;
  std::string _ls1ConfigFileName;
};

CPPUNIT_TEST_SUITE_REGISTRATION(LS1MoleculeIteratorTest);