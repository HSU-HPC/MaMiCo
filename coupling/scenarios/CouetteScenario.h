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

     _preciceAdapter = new coupling::solvers::PreciceAdapter<dim>();
    if (_preciceAdapter != NULL)
      std::cout << "MaMiCo: CouetteScenario::run(): Macro solver not NULL on rank " << _rank << std::endl;
    coupling::indexing::IndexingService<dim>::getInstance().init(_mdConfig, _mamicoConfig, _preciceAdapter, _rank);
    _multiMDCellService = new coupling::services::MultiMDCellService<simplemd::LinkedCell, dim>(_instanceHandling->getMDSolverInterface(), _preciceAdapter, _mdConfig, 
				_mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, dim>::getInstance().setMacroscopicSolverInterface(_preciceAdapter);
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    _multiMDCellService->computeAndStoreTemperature(_scenarioConfig.temp);
    allocateCFDToMDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(), *_preciceAdapter);
    allocateMDToCFDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(), *_preciceAdapter);
    double precice_dt = _preciceAdapter->initialize(
        _mdConfig.getDomainConfiguration().getGlobalDomainOffset(), _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
        _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(), _buf.globalCellIndices4MDToCFDBuffer, _buf.MDToCFDBuffer.size());
    std::cout << "MaMiCo: precice_dt=" << precice_dt << std::endl;
    std::cout << "MaMiCo: CouetteScenario::run() finish initialization " << std::endl;

    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      if (_preciceAdapter->isReadDataAvailable()) {
        _preciceAdapter->readData();
        const coupling::IndexConversion<3>& indexConversion{_multiMDCellService->getIndexConversion()};
        const unsigned int size = _buf.CFDToMDBuffer.size();
        const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
        const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
        for (unsigned int i = 0; i < size; i++) {
          const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf.globalCellIndices4CFDToMDBuffer[i]));
          tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
          for (unsigned int d = 0; d < 3; d++) {
            cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
          }
          double mass = macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
          tarch::la::Vector<3, double> momentum(mass * _preciceAdapter->getVelocity(cellMidPoint));
          _buf.CFDToMDBuffer[i]->setMicroscopicMass(mass);
          _buf.CFDToMDBuffer[i]->setMicroscopicMomentum(momentum);
        }
        _multiMDCellService->sendFromMacro2MD(_buf.CFDToMDBuffer, _buf.globalCellIndices4CFDToMDBuffer);
      }
      _instanceHandling->simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter, *_multiMDCellService);
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      _mdStepCounter += _mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
      double dt = _mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * _mdConfig.getSimulationConfiguration().getDt();
      if (_preciceAdapter->isWriteDataRequired(dt)) {
        _preciceAdapter->setMDBoundaryValues(_buf.MDToCFDBuffer, _buf.globalCellIndices4MDToCFDBuffer);
        _preciceAdapter->writeData();
      }
      precice_dt=_preciceAdapter->advance(dt);
      std::cout << "MaMiCo: precice_dt=" << precice_dt << std::endl;
      cycle++;
    }
    
    deleteBuffer(_buf.CFDToMDBuffer);
    if (_buf.globalCellIndices4CFDToMDBuffer != NULL) {
      delete[] _buf.globalCellIndices4CFDToMDBuffer;
      _buf.globalCellIndices4CFDToMDBuffer = NULL;
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
  void allocateCFDToMDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    const tarch::la::Vector<3, unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<3, unsigned int>(2));
    const unsigned int num = cells[0] * cells[1] * cells[2];
    deleteBuffer(_buf.CFDToMDBuffer);
    unsigned int numCells = 0;
    for (unsigned int i = 0; i < num; i++) {
      if (macroSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = macroSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          numCells++;
        }
      }
    }
    std::cout << "MaMiCo: CouetteScenario::allocateCFDToMDBuffer: numCells=" << numCells << std::endl;
    unsigned int* indices = new unsigned int[numCells];
    for (unsigned int i = 0; i < num; i++) {
      if (macroSolverInterface.sendMacroscopicQuantityToMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = macroSolverInterface.getSourceRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf.CFDToMDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          indices[_buf.CFDToMDBuffer.size() - 1] = i;
        }
      }
    }
    _buf.globalCellIndices4CFDToMDBuffer = indices;
  }

  void allocateMDToCFDBuffer(const coupling::IndexConversion<3>& indexConversion, coupling::interface::MacroscopicSolverInterface<3>& macroSolverInterface) {
    const tarch::la::Vector<3, unsigned int> cells(indexConversion.getGlobalNumberMacroscopicCells() + tarch::la::Vector<3, unsigned int>(2));
    const unsigned int num = cells[0] * cells[1] * cells[2];
    deleteBuffer(_buf.MDToCFDBuffer);
    unsigned int numCells = 0;
    for (unsigned int i = 0; i < num; i++) {
      if (macroSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = macroSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          numCells++;
        }
      }
    }
    std::cout << "MaMiCo: CouetteScenario::allocateMDToCFDBuffer: numCells=" << numCells << std::endl;
    unsigned int* indices = new unsigned int[numCells];
    for (unsigned int i = 0; i < num; i++) {
      if (macroSolverInterface.receiveMacroscopicQuantityFromMDSolver(indexConversion.getGlobalVectorCellIndex(i))) {
        std::vector<unsigned int> ranks = macroSolverInterface.getTargetRanks(indexConversion.getGlobalVectorCellIndex(i));
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf.MDToCFDBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          indices[_buf.MDToCFDBuffer.size() - 1] = i;
        }
      }
    }
    _buf.globalCellIndices4MDToCFDBuffer = indices;
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
    std::vector<coupling::datastructures::MacroscopicCell<3>*> CFDToMDBuffer;
    unsigned int* globalCellIndices4CFDToMDBuffer;
    std::vector<coupling::datastructures::MacroscopicCell<3>*> MDToCFDBuffer;
    unsigned int* globalCellIndices4MDToCFDBuffer;
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
