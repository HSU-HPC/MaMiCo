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
    init();
    for (int cycle = 0; cycle < _scenarioConfig.couplingCycles; cycle++) {
      advanceMacro(cycle);
      advanceMicro(cycle);
    }
    shutdown();
  }

private:
  void init() {
    const unsigned int dim = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    std::string xmlConfigurationFilename("config.xml");
    tarch::configuration::ParseConfiguration::parseConfiguration<simplemd::configurations::MolecularDynamicsConfiguration>(xmlConfigurationFilename,
                                                                                                                           "molecular-dynamics", _mdConfig);
    tarch::configuration::ParseConfiguration::parseConfiguration<coupling::configurations::MaMiCoConfiguration<dim>>(xmlConfigurationFilename, "mamico",
                                                                                                                     _mamicoConfig);
    _scenarioConfig = coupling::configurations::ScenarioConfig::parseConfiguration(xmlConfigurationFilename);
    _macroSolver = new coupling::solvers::PreciceAdapter<dim>(_mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0],
                                                              _mdConfig.getSimulationConfiguration().getDt() *
                                                                  _mdConfig.getSimulationConfiguration().getNumberOfTimesteps(),
                                                              _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap());
    if (_macroSolver != NULL)
      std::cout << "CouetteScenario::init() Macro solver not NULL on rank " << _rank << std::endl;
    _multiMDService = new tarch::utils::MultiMDService<dim>(_mdConfig.getMPIConfiguration().getNumberOfProcesses(), _scenarioConfig.totalNumberMDSimulations);
    _instanceHandling = new coupling::InstanceHandling<MY_LINKEDCELL, 3>(_mdConfig, _mamicoConfig, *_multiMDService);
    _mdStepCounter = 0;
    if (_rank == 0)
      gettimeofday(&_tv.start, NULL);
    _instanceHandling->switchOffCoupling();
    _instanceHandling->equilibrate(_scenarioConfig.equSteps, _mdStepCounter);
    _instanceHandling->switchOnCoupling();
    _mdStepCounter += _scenarioConfig.equSteps;
    _instanceHandling->setMDSolverInterface();
    coupling::indexing::IndexingService<dim>::getInstance().init(_mdConfig, _mamicoConfig, _macroSolver, _rank);
    _multiMDCellService = new coupling::services::MultiMDCellService<MY_LINKEDCELL, dim>(_instanceHandling->getMDSolverInterface(), _macroSolver, _mdConfig,
                                                                                         _mamicoConfig, xmlConfigurationFilename.c_str(), *_multiMDService);
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, dim>::getInstance().setMacroscopicSolverInterface(_macroSolver);
    _instanceHandling->setMacroscopicCellServices(*_multiMDCellService);
    _multiMDCellService->computeAndStoreTemperature(_scenarioConfig.temp);
    allocateCFDToMDBuffer(_multiMDCellService->getMacroscopicCellService(0).getIndexConversion(), *_macroSolver);
    _macroSolver->setCouplingMesh(_mdConfig.getDomainConfiguration().getGlobalDomainOffset(),
                                  _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize()[0], _buf.globalCellIndices4MDToCFDBuffer, _buf.MDToCFDBuffer.size());
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      double runtime = (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      std::cout << "Initialization: " << (int)(runtime / 1000) << "ms" << std::endl;
      gettimeofday(&_tv.start_total, NULL);
    }
    std::cout << "Finish CouetteScenario::init() " << std::endl;
  }

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
    std::cout << "CouetteScenario::allocateCFDToMDBuffer: numCells" << numCells << std::endl;
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
    std::cout << "CouetteScenario::allocateMDToCFDBuffer: numCells" << numCells << std::endl;
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

  void advanceMacro(int cycle) {
    if (_macroSolver != NULL) {
      if (_rank == 0) {
        gettimeofday(&_tv.start, NULL);
      }
      _macroSolver->advance(_mdConfig.getSimulationConfiguration().getDt() * _mdConfig.getSimulationConfiguration().getNumberOfTimesteps());
      if (_rank == 0) {
        gettimeofday(&_tv.end, NULL);
        _tv.macro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
      }
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
        tarch::la::Vector<3, double> momentum(mass * _macroSolver->getVelocity(cellMidPoint));
        _buf.CFDToMDBuffer[i]->setMicroscopicMass(mass);
        _buf.CFDToMDBuffer[i]->setMicroscopicMomentum(momentum);
      }
    }
    _multiMDCellService->sendFromMacro2MD(_buf.CFDToMDBuffer, _buf.globalCellIndices4CFDToMDBuffer);
  }

  void advanceMicro(int cycle) {
    if (_rank == 0) {
      gettimeofday(&_tv.start, NULL);
    }
    _instanceHandling->simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps(), _mdStepCounter, *_multiMDCellService);
    _multiMDCellService->plotEveryMacroscopicTimestep(cycle);
    _mdStepCounter += _mdConfig.getSimulationConfiguration().getNumberOfTimesteps();
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      _tv.micro += (_tv.end.tv_sec - _tv.start.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start.tv_usec);
    }
  }

  void shutdown() {
    if (_rank == 0) {
      gettimeofday(&_tv.end, NULL);
      double time_total = (_tv.end.tv_sec - _tv.start_total.tv_sec) * 1000000 + (_tv.end.tv_usec - _tv.start_total.tv_usec);
      std::cout << "Finished all coupling cycles after " << time_total / 1000000 << " s" << std::endl;
      std::cout << "Time percentages Micro, Macro, Other: " << std::endl;
      std::cout << _tv.micro / time_total * 100 << ", " << _tv.macro / time_total * 100 << ",  " << (1 - (_tv.micro + _tv.macro) / time_total) * 100
                << std::endl;
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
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, 3>::getInstance().getMacroscopicSolverInterface();
    if (_multiMDService != NULL) {
      delete _multiMDService;
      _multiMDService = NULL;
    }
    if (macroSolverInterface != NULL) {
      delete macroSolverInterface;
      macroSolverInterface = NULL;
    }
    if (_macroSolver != NULL) {
      delete _macroSolver;
      _macroSolver = NULL;
    }
    if (_multiMDCellService != NULL) {
      delete _multiMDCellService;
      _multiMDCellService = NULL;
    }
    std::cout << "Finish CouetteScenario::shutdown() " << std::endl;
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

  struct TimingValues {
    timeval start_total;
    timeval start;
    timeval end;
    double micro{0};
    double macro{0};
  };

  int _rank{0};
  simplemd::configurations::MolecularDynamicsConfiguration _mdConfig;
  coupling::configurations::MaMiCoConfiguration<3> _mamicoConfig;
  coupling::configurations::ScenarioConfig _scenarioConfig;
  unsigned int _mdStepCounter{0};
  coupling::solvers::PreciceAdapter<3>* _macroSolver;
  coupling::InstanceHandling<MY_LINKEDCELL, 3>* _instanceHandling;
  tarch::utils::MultiMDService<3>* _multiMDService;
  coupling::services::MultiMDCellService<MY_LINKEDCELL, 3>* _multiMDCellService;
  CouplingBuffer _buf;
  TimingValues _tv;
};

#endif // _COUPLING_TESTS_PRECICETEST_H_
