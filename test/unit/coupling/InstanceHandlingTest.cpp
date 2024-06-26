#include "coupling/InstanceHandling.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/CouetteConfiguration.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "simplemd/LinkedCell.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

/** tests several properties of the ErrorEstimation.
 *  @author Vahid Jafari
 */
class InstanceHandlingTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(InstanceHandlingTest);
  CPPUNIT_TEST(test);
  CPPUNIT_TEST_SUITE_END();

public:
  void test() {

    const char* filename("../test/unit/coupling/couette.xml");
    // std::string filename("../test/unit/coupling/couette.xml");
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename, "molecular-dynamics",
                                                                                                                           _simpleMDConfig);
    _cfg = coupling::configurations::CouetteConfig::parseCouetteConfiguration(filename);

    if (!_simpleMDConfig.isValid())
      CPPUNIT_FAIL("ERROR InstanceHandlingTest: Invalid SimpleMD config!");

    _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), _cfg.totalNumberMDSimulations
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                          ,
                                                          MPI_COMM_WORLD
#endif
    );

    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<3>>(filename, "mamico", _mamicoConfig);
    if (!_mamicoConfig.isValid())
      CPPUNIT_FAIL("ERROR InstanceHandlingTest: Invalid MaMiCo config!");
    _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(_simpleMDConfig, _mamicoConfig, *_multiMDService);
    if (_instanceHandling == nullptr)
      CPPUNIT_FAIL("ERROR InstanceHandlingTest::initSolvers() : _instanceHandling == NULL!");

    //    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations());
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getLocalNumberOfMDSimulations());

    //    CPPUNIT_ASSERT(_instanceHandling->getMDSolverInterface().size() != _multiMDService->getTotalNumberOfMDSimulations());
    CPPUNIT_ASSERT(_instanceHandling->getMDSolverInterface().size() != _multiMDService->getLocalNumberOfMDSimulations());
    _instanceHandling->setMDSolverInterface();
    //    CPPUNIT_ASSERT(_instanceHandling->getMDSolverInterface().size() == _multiMDService->getTotalNumberOfMDSimulations());
    CPPUNIT_ASSERT(_instanceHandling->getMDSolverInterface().size() == _multiMDService->getLocalNumberOfMDSimulations());
    // equilibrate MD
    _instanceHandling->switchOffCoupling();
    unsigned int t = 5;
    unsigned int T = 0;
    _instanceHandling->equilibrate(t, T);
    _instanceHandling->simulateTimesteps(t, T, 0);
    _instanceHandling->simulateTimesteps(t, T, 1);
    _instanceHandling->simulateTimesteps(t, T);

    _instanceHandling->switchOnCoupling();
    //      unsigned int  totalNumberOfMDSimulations = _multiMDService->getTotalNumberOfMDSimulations();
    _instanceHandling->addSimulationBlock();
    //      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations()+1);
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getLocalNumberOfMDSimulations() + 1);
    _instanceHandling->rmSimulationBlock();
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getLocalNumberOfMDSimulations());
    _instanceHandling->addSimulationBlock();
    auto* mdSolverInterface =
        _instanceHandling->addMDSimulation(_multiMDService->getLocalNumberOfMDSimulations(), _multiMDService->getLocalNumberOfMDSimulations());
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getLocalNumberOfMDSimulations() + 1);
    CPPUNIT_ASSERT(mdSolverInterface != nullptr);
    _instanceHandling->rmMDSimulation(_instanceHandling->getSimpleMD().size() - 1);
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD()[_instanceHandling->getSimpleMD().size() - 1] == nullptr);
    _instanceHandling->rmSimulationBlock();
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getLocalNumberOfMDSimulations());

    _instanceHandling->writeCheckpoint("restart_checkpoint", 0);
    std::ifstream inputFile("restart_checkpoint_0_0.checkpoint");
    CPPUNIT_ASSERT(inputFile.is_open());
    std::cout << "restart_checkpoint_0_0.checkpoint read!" << std::endl;

    _instanceHandling->~InstanceHandling();
    CPPUNIT_ASSERT(!_instanceHandling->getSimpleMD().size());
    CPPUNIT_ASSERT(!_instanceHandling->getMDSolverInterface().size());
  }

private:
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
  coupling::configurations::CouetteConfig _cfg;
};

CPPUNIT_TEST_SUITE_REGISTRATION(InstanceHandlingTest);
