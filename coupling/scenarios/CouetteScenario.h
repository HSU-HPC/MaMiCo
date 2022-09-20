#ifndef _COUPLING_TESTS_PRECICETEST_H_
#define _COUPLING_TESTS_PRECICETEST_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/configurations/ScenarioConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/scenarios/Scenario.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/solvers/CouetteSolverInterface.h"
#include "coupling/solvers/PreciceAdapter.h"
#include "precice/SolverInterface.hpp"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include <mpi.h>
#include <random>
#include <sys/time.h>

class CouetteScenario : public Scenario {
public:
  CouetteScenario() : Scenario("CouetteScenario") {}
  virtual ~CouetteScenario() {}

  virtual void run() {
    const unsigned int dim = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    std::string xmlConfigurationFilename("config.xml");
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename,
                                                                                                                           "molecular-dynamics", _mdConfig);
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico",
                                                                                                                     _mamicoConfig);
    _scenarioConfig = coupling::configurations::ScenarioConfig::parseConfiguration(xmlConfigurationFilename);

    _multiMDService = new tarch::utils::MultiMDService<dim>(_mdConfig.getMPIConfiguration().getNumberOfProcesses(), _scenarioConfig.totalNumberMDSimulations);
    _instanceHandling = new coupling::InstanceHandling<simplemd::LinkedCell, 3>(_mdConfig, _mamicoConfig, *_multiMDService);
    _mdStepCounter = 0;
    _instanceHandling->switchOffCoupling();
    _instanceHandling->equilibrate(_scenarioConfig.equSteps, _mdStepCounter);
    _instanceHandling->switchOnCoupling();
    _mdStepCounter += _scenarioConfig.equSteps;
    _instanceHandling->setMDSolverInterface();

     _preciceAdapter = new coupling::solvers::PreciceAdapter<dim>(_mdConfig.getDomainConfiguration().getGlobalDomainOffset(), 
     _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(), _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
    if (_preciceAdapter != NULL)
      std::cout << "MaMiCo: CouetteScenario::run(): Macro solver not NULL on rank " << _rank << std::endl;
    coupling::indexing::IndexingService<dim>::getInstance().init(_mdConfig, _mamicoConfig, _preciceAdapter, _rank);
    _multiMDCellService = new coupling::services::MultiMDCellService<simplemd::LinkedCell, dim>(_instanceHandling->getMDSolverInterface(), _preciceAdapter, _mdConfig, 
				_mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, dim>::getInstance().setMacroscopicSolverInterface(_preciceAdapter);
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    _multiMDCellService->computeAndStoreTemperature(_scenarioConfig.temp);
    allocateM2mBuffer();
    allocatem2MBuffer();
    double precice_dt = _preciceAdapter->initialize(
        _buf.M2mCellGlobalIndices, _buf.M2mBuffer.size(),
        _buf.m2MCellGlobalIndices, _buf.m2MBuffer.size());
    std::cout << "MaMiCo: CouetteScenario::run() finish initialization " << std::endl;
    std::cout << "MaMiCo: precice_dt=" << precice_dt << std::endl;
    double mamico_dt = _mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * _mdConfig.getSimulationConfiguration().getDt();
    std::cout << "MaMiCo: mamico_dt=" << mamico_dt << std::endl;
    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      if (_preciceAdapter->isReadDataAvailable()) {
        _preciceAdapter->readData();
        const coupling::IndexConversion<3>& indexConversion{_multiMDCellService->getIndexConversion()};
        const unsigned int size = _buf.M2mBuffer.size();
        const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
        const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
        for (unsigned int i = 0; i < size; i++) {
          const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf.M2mCellGlobalIndices[i]));
          tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
          for (unsigned int d = 0; d < 3; d++) {
            cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
          }
          double mass = macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
          tarch::la::Vector<3, double> momentum(mass * _preciceAdapter->getVelocity(cellMidPoint));
          _buf.M2mBuffer[i]->setMicroscopicMass(mass);
          _buf.M2mBuffer[i]->setMicroscopicMomentum(momentum);
        }
        _multiMDCellService->sendFromMacro2MD(_buf.M2mBuffer, _buf.M2mCellGlobalIndices);
      }
      // should be like that
      //double dt = std::min(precice_dt, mamico_dt);
      std::cout << "MaMiCo: precice_dt=" << precice_dt << std::endl;
      std::cout << "MaMiCo: mamico_dt=" << mamico_dt << std::endl;
      _instanceHandling->simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter, *_multiMDCellService);
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      _mdStepCounter += _mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
      if (_preciceAdapter->isWriteDataRequired(mamico_dt)) {
        _preciceAdapter->setMDBoundaryValues(_buf.m2MBuffer, _buf.m2MCellGlobalIndices);
        _preciceAdapter->writeData();
      }
      precice_dt=_preciceAdapter->advance(mamico_dt);
      cycle++;
    }
    
    deleteBuffer(_buf.M2mBuffer);
    if (_buf.M2mCellGlobalIndices != NULL) {
      delete[] _buf.M2mCellGlobalIndices;
      _buf.M2mCellGlobalIndices = NULL;
    }
    if (_instanceHandling != nullptr) {
      delete _instanceHandling;
    }
    coupling::interface::MacroscopicSolverInterface<3>* macroSolverInterface =
        coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, 3>::getInstance().getMacroscopicSolverInterface();
    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (macroSolverInterface != NULL) {
      delete macroSolverInterface;
      macroSolverInterface = NULL;
    }
    if (_preciceAdapter != NULL) {
      delete _preciceAdapter;
      _preciceAdapter = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
    std::cout << "MaMiCo: Finish CouetteScenario::shutdown() " << std::endl;
  }
  

private:
  void allocateM2mBuffer() {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> M2mBufferCellGlobalIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_preciceAdapter->sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _preciceAdapter->getSourceRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf.M2mBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          M2mBufferCellGlobalIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    _buf.M2mCellGlobalIndices = M2mBufferCellGlobalIndices.data();
    std::cout << "MaMiCo: CouetteScenario::allocateM2mBuffer: numCells=" << _buf.M2mBuffer.size() << std::endl;
  }

  void allocatem2MBuffer() {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> m2MBufferCellGlobalIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_preciceAdapter->receiveMacroscopicQuantityFromMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _preciceAdapter->getTargetRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf.m2MBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          m2MBufferCellGlobalIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    _buf.m2MCellGlobalIndices = m2MBufferCellGlobalIndices.data();
    std::cout << "MaMiCo: CouetteScenario::allocatem2MBuffer: numCells=" << _buf.m2MBuffer.size() << std::endl;
  }

  void deleteBuffer(std::vector<coupling::datastructures::MacroscopicCell<3>*>& buffer) const {
    for (unsigned int i = 0; i < buffer.size(); i++) {
      if (buffer[i] != NULL) {
        delete buffer[i];
        buffer[i] = NULL;
      }
    }
    buffer.clear();
  }

  struct CouplingBuffer {
    std::vector<coupling::datastructures::MacroscopicCell<3>*> M2mBuffer;
    unsigned int* M2mCellGlobalIndices;
    std::vector<coupling::datastructures::MacroscopicCell<3>*> m2MBuffer;
    unsigned int* m2MCellGlobalIndices;
  };

  int _rank{0};
  simplemd::configurations::MolecularDynamicsConfiguration _mdConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  coupling::configurations::ScenarioConfig _scenarioConfig;
  unsigned int _mdStepCounter{0};
  coupling::solvers::PreciceAdapter<3>* _preciceAdapter;
  coupling::InstanceHandling<simplemd::LinkedCell, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<simplemd::LinkedCell, 3>* _multiMDCellService;
  CouplingBuffer _buf;
};

#endif // _COUPLING_TESTS_PRECICETEST_H_
