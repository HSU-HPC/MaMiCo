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
#include "coupling/scenarios/PreciceAdapter.h"
#include "precice/SolverInterface.hpp"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/utils/MultiMDService.h"
#include <mpi.h>
#include <random>
#include <sys/time.h>

class CouetteScenario : public Scenario {
public:
  CouetteScenario() : Scenario("CouetteScenario"), _logger(tarch::logging::Logger("CouetteScenario")) {}
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

    _preciceAdapter = new PreciceAdapter<dim>(_mdConfig.getDomainConfiguration().getGlobalDomainOffset(), _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize());
    tarch::la::Vector<3, int> globalNumberMacroscopicCells;
    for (unsigned int d = 0; d < 3; d++) {
      globalNumberMacroscopicCells[d] = floor(_mdConfig.getDomainConfiguration().getGlobalDomainSize()[d] / _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[d] + 0.5);
    }
    _macroscopicSolverInterface = new MyMacroscopicSolverInterface<dim>(globalNumberMacroscopicCells, (int)_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap(),  _rank);
    coupling::indexing::IndexingService<dim>::getInstance().init(_mdConfig, _mamicoConfig, _macroscopicSolverInterface, _rank);
    _multiMDCellService = new coupling::services::MultiMDCellService<simplemd::LinkedCell, dim>(_instanceHandling->getMDSolverInterface(), _macroscopicSolverInterface, _mdConfig, 
				_mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, dim>::getInstance().setMacroscopicSolverInterface(_macroscopicSolverInterface);
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    _multiMDCellService->computeAndStoreTemperature(_scenarioConfig.temp);
    allocateMacro2MicroBuffer();
    allocateMicro2MacroBuffer();
    double precice_dt = _preciceAdapter->initialize(
        _buf._macro2MicroCellGlobalIndices, _buf._macro2MicroBuffer.size()
#ifdef TWO_WAY 
        , _buf._micro2MacroCellGlobalIndices, _buf._micro2MacroBuffer.size());
#else
        );
#endif
    _logger.info("rank {} finished initialization", _rank);
    _logger.info("rank {} precice_dt={}", _rank, precice_dt);
    double mamico_dt = _mdConfig.getSimulationConfiguration().getNumberOfTimesteps() * _mdConfig.getSimulationConfiguration().getDt();
    _logger.info("rank {} precice_dt={}", _rank, mamico_dt);
    int cycle = 0;
    while (_preciceAdapter->isCouplingOngoing()) {
      _preciceAdapter->readData();
      // _logger.info("rank {} data read from precice", _rank);
      const coupling::IndexConversion<3>& indexConversion{_multiMDCellService->getIndexConversion()};
      const unsigned int size = _buf._macro2MicroBuffer.size();
      const tarch::la::Vector<3, double> domainOffset(indexConversion.getGlobalMDDomainOffset());
      const tarch::la::Vector<3, double> macroscopicCellSize(indexConversion.getMacroscopicCellSize());
      for (unsigned int i = 0; i < size; i++) {
        const tarch::la::Vector<3, unsigned int> globalIndex(indexConversion.getGlobalVectorCellIndex(_buf._macro2MicroCellGlobalIndices[i]));
        tarch::la::Vector<3, double> cellMidPoint(domainOffset - 0.5 * macroscopicCellSize);
        for (unsigned int d = 0; d < 3; d++) {
          cellMidPoint[d] = cellMidPoint[d] + ((double)globalIndex[d]) * macroscopicCellSize[d];
        }
        double mass = macroscopicCellSize[0] * macroscopicCellSize[1] * macroscopicCellSize[2];
        tarch::la::Vector<3, double> momentum(mass * _preciceAdapter->getVelocity(cellMidPoint));
        _buf._macro2MicroBuffer[i]->setMicroscopicMass(mass);
        _buf._macro2MicroBuffer[i]->setMicroscopicMomentum(momentum);
      } 
      _multiMDCellService->sendFromMacro2MD(_buf._macro2MicroBuffer, _buf._macro2MicroCellGlobalIndices);
      // should be like that
      //double dt = std::min(precice_dt, mamico_dt);
      // _logger.info("rank {} precice_dt={}", _rank, precice_dt);
      // _logger.info("rank {} precice_dt={}", _rank, mamico_dt);
      _instanceHandling->simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter, *_multiMDCellService);
      _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
      _mdStepCounter += _mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
      _multiMDCellService->sendFromMD2Macro(_buf._micro2MacroBuffer, _buf._micro2MacroCellGlobalIndices);
#ifdef TWO_WAY 
      _preciceAdapter->setMDBoundaryValues(_buf._micro2MacroBuffer, _buf._micro2MacroCellGlobalIndices);
      _preciceAdapter->writeData();
#endif
      precice_dt=_preciceAdapter->advance(mamico_dt);
      cycle++;
      if (_buf._micro2MacroBuffer.size() != 0 && _scenarioConfig.csvEveryTimestep >= 1 && cycle % _scenarioConfig.csvEveryTimestep == 0)
        write2CSV(_buf._micro2MacroBuffer, _buf._micro2MacroCellGlobalIndices, cycle, globalNumberMacroscopicCells, (int)_mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
    }
    
    deleteBuffer(_buf._macro2MicroBuffer);
    if (_buf._macro2MicroCellGlobalIndices != NULL) {
      delete[] _buf._macro2MicroCellGlobalIndices;
      _buf._macro2MicroCellGlobalIndices = NULL;
    }
    deleteBuffer(_buf._micro2MacroBuffer);
    if (_buf._micro2MacroCellGlobalIndices != NULL) {
      delete[] _buf._micro2MacroCellGlobalIndices;
      _buf._micro2MacroCellGlobalIndices = NULL;
    }
    if (_instanceHandling != nullptr) {
      delete _instanceHandling;
    }
    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (_macroscopicSolverInterface != NULL) {
      delete _macroscopicSolverInterface;
      _macroscopicSolverInterface = NULL;
    }
    if (_preciceAdapter != NULL) {
      delete _preciceAdapter;
      _preciceAdapter = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
    _logger.info("rank {} shuttdown", _rank);
  }
  

private:
  void allocateMacro2MicroBuffer() {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> macro2MicroBufferCellGlobalIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_macroscopicSolverInterface->sendMacroscopicQuantityToMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getSourceRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf._macro2MicroBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          macro2MicroBufferCellGlobalIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    _buf._macro2MicroCellGlobalIndices = new unsigned int[_buf._macro2MicroBuffer.size()];
    std::copy(macro2MicroBufferCellGlobalIndices.begin(), macro2MicroBufferCellGlobalIndices.end(), _buf._macro2MicroCellGlobalIndices);
    _logger.debug("rank {} macro2Micro buffer numCells {}", _rank, _buf._macro2MicroBuffer.size());
  }

  void allocateMicro2MacroBuffer() {
    using namespace coupling::indexing;
    const unsigned int dim = 3;
    std::vector<unsigned int> micro2MacroBufferCellGlobalIndices;
    for (CellIndex<dim, IndexTrait::vector> cellIndex : CellIndex<dim, IndexTrait::vector>()) {
      tarch::la::Vector<3, unsigned int> cellVectorIndex = static_cast<tarch::la::Vector<dim, unsigned int>>(cellIndex.get());
      if (_macroscopicSolverInterface->receiveMacroscopicQuantityFromMDSolver(cellVectorIndex)) {
        std::vector<unsigned int> ranks = _macroscopicSolverInterface->getTargetRanks(cellVectorIndex);
        bool containsThisRank = false;
        for (unsigned int k = 0; k < ranks.size(); k++) {
          containsThisRank = containsThisRank || (ranks[k] == (unsigned int)_rank);
        }
        if (containsThisRank) {
          _buf._micro2MacroBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
          micro2MacroBufferCellGlobalIndices.push_back(convertToScalar<dim>(cellIndex));
        }
      }
    }
    _buf._micro2MacroCellGlobalIndices = new unsigned int[_buf._micro2MacroBuffer.size()];
    std::copy(micro2MacroBufferCellGlobalIndices.begin(), micro2MacroBufferCellGlobalIndices.end(), _buf._micro2MacroCellGlobalIndices);
    _logger.debug("rank {} micro2Macro buffer numCells {}", _rank, _buf._micro2MacroBuffer.size());
  }

  void write2CSV(std::vector<coupling::datastructures::MacroscopicCell<3>*>& micro2MacroBuffer, const unsigned int* const micro2MacroCellGlobalIndices,
                 int couplingCycle, const tarch::la::Vector<3, int> globalNumberMacroscopicCells, const int overlap) const {
    std::stringstream ss;
    ss << "results_" << _rank << "_" << couplingCycle << ".csv";
    std::ofstream file(ss.str().c_str());
    if (!file.is_open()) {
      exit(EXIT_FAILURE);
    }
    tarch::la::Vector<3, int> cellIndex;
    unsigned int count = 0;
    for (unsigned int i = 0; i < micro2MacroBuffer.size(); i++) {
      cellIndex = coupling::indexing::convertToVector<3>({micro2MacroCellGlobalIndices[i]});
      bool isOuter = false;
      for (unsigned int currentDim = 0; currentDim < 3; currentDim++) {
        isOuter |= (int)cellIndex[currentDim] == overlap;
        isOuter |= (int)cellIndex[currentDim] == globalNumberMacroscopicCells[currentDim] + 1 - overlap;
      }
      if (!isOuter) {
        tarch::la::Vector<3, double> vel(micro2MacroBuffer[i]->getMacroscopicMomentum());
        if (micro2MacroBuffer[i]->getMacroscopicMass() != 0.0) {
          vel = (1.0 / micro2MacroBuffer[i]->getMacroscopicMass()) * vel;
        }
        file << cellIndex[0] << " ; " << cellIndex[1] << " ; " << cellIndex[2] << " ; " << vel[0] << " ; " << vel[1] << " ; " << vel[2] << " ; "
            << micro2MacroBuffer[i]->getMacroscopicMass() << ";";
        file << std::endl;
        count++;
      }
    }
    file.close();
    _logger.debug("Number of cells written = {}", count);
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
    std::vector<coupling::datastructures::MacroscopicCell<3>*> _macro2MicroBuffer;
    unsigned int* _macro2MicroCellGlobalIndices;
    std::vector<coupling::datastructures::MacroscopicCell<3>*> _micro2MacroBuffer;
    unsigned int* _micro2MacroCellGlobalIndices;
  };

  template <unsigned int dim> class MyMacroscopicSolverInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
  public:
    MyMacroscopicSolverInterface(const tarch::la::Vector<3, int> globalNumberMacroscopicCells, const int overlap, const int rank) : _globalNumberMacroscopicCells(globalNumberMacroscopicCells), _overlap(overlap), _rank(rank) {}

    bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
      bool rcv = true;
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
        rcv &= (int)globalCellIndex[currentDim] >= 1 + (_overlap - 1) ;
        rcv &= (int)globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - (_overlap - 1);
      }
      return rcv;
    }

    bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) override {
      bool isGhostCell = false;
      bool isInner = true;
      for (unsigned int currentDim = 0; currentDim < dim; currentDim++) {
        isGhostCell |= (int)globalCellIndex[currentDim] > _globalNumberMacroscopicCells[currentDim];
        isGhostCell |= (int)globalCellIndex[currentDim] < 1;
        isInner &= (int)globalCellIndex[currentDim] >= 1 + _overlap;
        isInner &= (int)globalCellIndex[currentDim] < _globalNumberMacroscopicCells[currentDim] + 1 - _overlap;
      }
      return (!isGhostCell) && (!isInner);
    }

    std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
      return {0};
    }

    std::vector<unsigned int> getSourceRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
      return {0};
    }

    std::vector<unsigned int> getTargetRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) override {
      return {0};
    }

  private:
    const tarch::la::Vector<3, int> _globalNumberMacroscopicCells;
    const int _overlap;
    const int _rank;
  };

  tarch::logging::Logger _logger;
  int _rank{0};
  simplemd::configurations::MolecularDynamicsConfiguration _mdConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  coupling::configurations::ScenarioConfig _scenarioConfig;
  unsigned int _mdStepCounter{0};
  PreciceAdapter<3>* _preciceAdapter;
  coupling::interface::MacroscopicSolverInterface<3>* _macroscopicSolverInterface;
  coupling::InstanceHandling<simplemd::LinkedCell, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<simplemd::LinkedCell, 3>* _multiMDCellService;
  CouplingBuffer _buf;
};

#endif // _COUPLING_TESTS_PRECICETEST_H_