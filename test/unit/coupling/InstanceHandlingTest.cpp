#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "simplemd/LinkedCell.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/InstanceHandling.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"


/** tests several properties of the ErrorEstimation. 
 *  @author Vahid Jafarimann
 */
class InstanceHandlingTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(InstanceHandlingTest);
  CPPUNIT_TEST(test);
  CPPUNIT_TEST_SUITE_END();

public:
    
    void test() {
        
        std::string filename("couette.xml");
        tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(filename, "molecular-dynamics",_simpleMDConfig);
        
    if (!_simpleMDConfig.isValid()) {
      std::cout << "ERROR CouetteScenario: Invalid SimpleMD config!" << std::endl;
      exit(EXIT_FAILURE);
    }
        
        _multiMDService = new tarch::utils::MultiMDService<3>(_simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(), 2
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
,1
#endif
    );
        _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(_simpleMDConfig, _mamicoConfig, *_multiMDService);
        if (_instanceHandling == nullptr) {
      std::cout << "ERROR CouetteScenario::initSolvers() : _instanceHandling == NULL!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations());
    // equilibrate MD
      unsigned int t = 1;
      _instanceHandling->switchOffCoupling();
      _instanceHandling->equilibrate(2, 0);
      _instanceHandling->switchOnCoupling();
//      _instanceHandling->simulateTimesteps(t, t);
      _instanceHandling->writeCheckpoint("check",10);
      std::ifstream inputFile("check_10_0.checkpoint");
      CPPUNIT_ASSERT(inputFile.is_open());
      unsigned int  totalNumberOfMDSimulations = _multiMDService->getTotalNumberOfMDSimulations();
      _instanceHandling->addSimulationBlock();
      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations()+1);
      _instanceHandling->rmSimulationBlock();
      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations());
      _instanceHandling->addSimulationBlock();
      auto* mdSolverInterface = _instanceHandling->addMDSimulation(totalNumberOfMDSimulations, totalNumberOfMDSimulations);
      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations()+1);
      CPPUNIT_ASSERT(mdSolverInterface != nullptr);
      _instanceHandling->rmMDSimulation(_instanceHandling->getSimpleMD().size()-1);
      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD()[_instanceHandling->getSimpleMD().size()-1] == nullptr);
      _instanceHandling->rmSimulationBlock();
      CPPUNIT_ASSERT(_instanceHandling->getSimpleMD().size() == _multiMDService->getTotalNumberOfMDSimulations());
      
      
      
      
      
      
      
      
      _instanceHandling->~InstanceHandling();   
    
  }
 

private:
    
    coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
    tarch::utils::MultiMDService<3>* _multiMDService;
    coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
    simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
 
 
};

CPPUNIT_TEST_SUITE_REGISTRATION(InstanceHandlingTest);
