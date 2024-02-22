// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_

namespace coupling {
template <class LinkedCell, unsigned int dim> class InstanceHandling;
namespace services {
template <class LinkedCell, unsigned int dim> class MultiMDCellService;
}
} // namespace coupling

#include "coupling/InstanceHandling.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/services/CouplingCellService.h"
#include "coupling/services/CouplingCellServiceDummy.h"

/** service to couple various MD simulations to a single instance of a
 * macroscopic solver. We currently consider only one layout: global ranks are
 * linearly split into blocks of size numberProcesses. Each block may run
 * several MD simulations. Each of these simulations has one mdSolverInterface
 * (handed over in the constructor as vector).
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::services::MultiMDCellService {
public:
  MultiMDCellService(std::vector<coupling::interface::MDSolverInterface<LinkedCell, dim>*> mdSolverInterfaces, // MD solver interfaces for each MD simulation
                                                                                                               // (which uses this rank)
                     coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface,
                     simplemd::configurations::MolecularDynamicsConfiguration& mdConfiguration,
                     coupling::configurations::MaMiCoConfiguration<dim>& mamicoConfiguration, const char* filterPipelineConfiguration,
                     tarch::utils::MultiMDService<dim>& multiMDService, int tws = 0)
      : _multiMDService(multiMDService), _tws(tws), _intNumberProcesses(computeScalarNumberProcesses()), _mdConfiguration(mdConfiguration),
        _mamicoConfiguration(mamicoConfiguration), _filterPipelineConfiguration(filterPipelineConfiguration),
        _macroscopicSolverInterface(macroscopicSolverInterface), _postMultiInstanceFilterPipeline(nullptr) {
    _topologyOffset = computeTopologyOffset();
    _localNumberMDSimulations = multiMDService.getLocalNumberOfMDSimulations();
    _totalNumberMDSimulations = multiMDService.getTotalNumberOfMDSimulations();

    _blockOffset = _localNumberMDSimulations * _topologyOffset / _intNumberProcesses;

    _listActiveMDSimulations = std::vector<bool>(_totalNumberMDSimulations, true);
    _nextFreeBlock = _multiMDService.getNumberLocalComms() - 1;
    _warmupPhase = std::vector<int>(_totalNumberMDSimulations, 0);

    // If we allow for zero MD instances on ranks, this initializion process
    // would segfault...
    if (mdSolverInterfaces.size() == 0) {
      std::cout << "ERROR: Zero MD instances on rank " << _multiMDService.getGlobalRank() << ". Maybe you forgot to adjust number-md-simulations?" << std::endl;
      exit(EXIT_FAILURE);
    }

    // determine globally unique IDs for each coupling cell service.
    // Assumptions:
    // - the global parallel topology is subdivided into equally sized blocks of
    // ranks
    // - on each block of ranks, the SAME number of MD simulations is executed.
    // This yields that the variable _topologyOffset together with the local ID
    // of an MD simulation yields a unique global ID
    // - an MD simulation executes on one full block of ranks
    _couplingCellServices = new coupling::services::CouplingCellService<dim>*[_totalNumberMDSimulations];
    if (_couplingCellServices == NULL) {
      std::cout << "ERROR "
                   "coupling::services::MultiMDCellService::MultiMDCellService("
                   "...): _couplingCellServices==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // allocate all coupling cell services for macro-only-solver BEFORE
    // topology offset
    // -> we have localNumberMDSimulations that run per block of ranks on
    // intNumberProcesses
    //    -> this yields
    // localNumberMDSimulations*topologyOffset/intNumberProcesses MD simulations
    // before the actual block of ranks
    for (unsigned int i = 0; i < _blockOffset; i++) {
      _couplingCellServices[i] = createCouplingCellServiceDummy(i, macroscopicSolverInterface, mdConfiguration, _multiMDService, mamicoConfiguration,
                                                                (i / _localNumberMDSimulations) * _intNumberProcesses);
      if (_couplingCellServices[i] == NULL) {
        std::cout << "ERROR "
                     "coupling::services::MultiMDCellService::MultiMDCellServic"
                     "e(...): _couplingCellServices["
                  << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    // allocate all coupling cell services for macro- and micro-interactions
    for (unsigned int i = _blockOffset; i < _localNumberMDSimulations + _blockOffset; i++) {
      _couplingCellServices[i] = new coupling::services::CouplingCellServiceImpl<LinkedCell, dim>(
          i, mdSolverInterfaces[i - _blockOffset], macroscopicSolverInterface, mdConfiguration.getMPIConfiguration().getNumberOfProcesses(),
          _multiMDService.getGlobalRank(), mamicoConfiguration.getParticleInsertionConfiguration(), mamicoConfiguration.getMomentumInsertionConfiguration(),
          mamicoConfiguration.getBoundaryForceConfiguration(), mamicoConfiguration.getTransferStrategyConfiguration(),
          mamicoConfiguration.getParallelTopologyConfiguration(), mamicoConfiguration.getThermostatConfiguration(),
          mdConfiguration.getSimulationConfiguration().getNumberOfTimesteps(), mamicoConfiguration.getCouplingCellConfiguration(),
          _filterPipelineConfiguration.c_str(), multiMDService, _topologyOffset, _tws);
      if (_couplingCellServices[i] == NULL) {
        std::cout << "ERROR "
                     "coupling::services::MultiMDCellService::MultiMDCellServic"
                     "e(...): _couplingCellServices["
                  << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    // allocate all coupling cell services for macro-only-solver AFTER
    // topology offset
    for (unsigned int i = _blockOffset + _localNumberMDSimulations; i < _totalNumberMDSimulations; i++) {
      _couplingCellServices[i] = createCouplingCellServiceDummy(i, macroscopicSolverInterface, mdConfiguration, _multiMDService, mamicoConfiguration,
                                                                (i / _localNumberMDSimulations) * _intNumberProcesses);
      if (_couplingCellServices[i] == NULL) {
        std::cout << "ERROR "
                     "coupling::services::MultiMDCellService::MultiMDCellServic"
                     "e(...): _couplingCellServices["
                  << i << "]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    _mdConfiguration.getDomainConfigurationNonConst().setInitFromCheckpoint(true);
    std::stringstream filestem;
    filestem << "restart_checkpoint_" << (_multiMDService.getGlobalRank()) / computeScalarNumberProcesses() << "_0";
    _mdConfiguration.getDomainConfigurationNonConst().setCheckpointFilestem(filestem.str());
    _mdConfiguration.getDomainConfigurationNonConst().setInitFromSequentialCheckpoint(false);
  }

  ~MultiMDCellService() {
    for (unsigned int i = 0; i < _totalNumberMDSimulations; i++) {
      if (_couplingCellServices[i] != nullptr) {
        delete _couplingCellServices[i];
        _couplingCellServices[i] = nullptr;
      }
    }
    if (_couplingCellServices != NULL) {
      delete[] _couplingCellServices;
      _couplingCellServices = NULL;
    }

    // TODO: fix free bug/possible memory leak here
    for (auto coupling_cell : _couplingCells) {
      delete coupling_cell;
    }
    _couplingCells.clear();
    // delete [] _couplingCells.data();
    if (_postMultiInstanceFilterPipeline != nullptr)
      delete _postMultiInstanceFilterPipeline;
  }

  /** get access to coupling cell services. This is required potentially by
   * each MD simulation to incorporate cell-local thermostat, mass and momentum
   * transfer */
  coupling::services::CouplingCellService<dim>& getCouplingCellService(unsigned int localIndex) {
    if (localIndex < _localNumberMDSimulations) {
      return *(_couplingCellServices[_topologyOffset * _localNumberMDSimulations / _intNumberProcesses + localIndex]);
    } else {
      std::cout << "ERROR MultiMDCellService::getCouplingCellService(localIndex): "
                   "localIndex >_localNumberMDSimulations-1!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** Forward Plotting to coupling cell services if not NULL! */
  void plotEveryMacroscopicTimestepforCouplingCellService(unsigned int localIndex, int cycle) {
    if (_couplingCellServices[_topologyOffset * _localNumberMDSimulations / _intNumberProcesses + localIndex] == nullptr)
      return;
    _couplingCellServices[_topologyOffset * _localNumberMDSimulations / _intNumberProcesses + localIndex]->plotEveryMacroscopicTimestep(cycle);
  }

  void plotEveryMacroscopicTimestep(int cycle) {
    for (unsigned int i = _blockOffset; i < _blockOffset + _localNumberMDSimulations; ++i) {
      if (_couplingCellServices[i] == nullptr)
        continue;
      _couplingCellServices[i]->plotEveryMacroscopicTimestep(cycle);
    }
  }

  void computeAndStoreTemperature(double temp) {
    for (unsigned int i = _blockOffset; i < _blockOffset + _localNumberMDSimulations; ++i) {
      if (_couplingCellServices[i] == nullptr)
        continue;
      _couplingCellServices[i]->computeAndStoreTemperature(temp);
    }
  }

  /** send data from macroscopic solver to all MD simulations */
  void sendFromMacro2MD(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver, const I00* const indices) {
    // just send information to all MD instances. This is currently
    // sequentialized inside CouplingCellService/SendRecvBuffer
    for (unsigned int i = 0; i < _totalNumberMDSimulations; i++) {
      // std::cout << "Rank " <<
      //_couplingCellServices[i]->getIndexConversion().getThisRank() << ":
      // Send from macro to MD for Simulation no. " << i << std::endl;
      if (_couplingCellServices[i] != nullptr) {
        _couplingCellServices[i]->sendFromMacro2MD(couplingCellsFromMacroscopicSolver, indices);
      }
    }
  }

  /** broadcasts data from macroscopic solver to all MD simulations */
  void bcastFromMacro2MD(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver, const I00* const indices) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_NO)
    // Fall back on sequential operation when MPI is not available (avoids redundant implementation)
    sendFromMacro2MD(couplingCellsFromMacroscopicSolver, indices);
    return;
#else

    std::vector<coupling::sendrecv::DataExchangeFromMacro2MD<dim>*> allDEs(_totalNumberMDSimulations);
    std::vector<std::vector<coupling::datastructures::CouplingCell<dim>*>> allCouplingCellsFromMamico(_totalNumberMDSimulations);
    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (nullptr == _couplingCellServices[i])
        continue;
      allDEs[i] = new coupling::sendrecv::DataExchangeFromMacro2MD<dim>(_macroscopicSolverInterface, _couplingCellServices[i]->getID());
      if (auto* v = dynamic_cast<CouplingCellServiceImpl<LinkedCell, dim>*>(_couplingCellServices[i])) {
        allCouplingCellsFromMamico[i] = v->getCouplingCells().getCouplingCells();
      }
    }

    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (nullptr == _couplingCellServices[i])
        continue;
      _couplingCellServices[i]->sendFromMacro2MDPreProcess();
    }

    coupling::sendrecv::FromMacro2MD<coupling::datastructures::CouplingCell<dim>, dim> fromMacro2MD;
    fromMacro2MD.bcastFromMacro2MD(allDEs, couplingCellsFromMacroscopicSolver, indices, allCouplingCellsFromMamico);

    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (nullptr == _couplingCellServices[i])
        continue;
      _couplingCellServices[i]->sendFromMacro2MDPostProcess();
    }
    for (coupling::sendrecv::DataExchangeFromMacro2MD<dim>*& de : allDEs) {
      delete de;
      de = nullptr;
    }
#endif
  }

  /** Creates the sum over all instances' coupling cells, in order to reduce the amount of communication needed.
   */
  void sumUpCouplingCellsFromMamico() {
    /*
     *  On the first coupling step,  initialize the reduced couplingCell vector according to `size`
     */
    if (!_sumCouplingCells.empty()) {
      for (unsigned int i = 0; i < _sumCouplingCells.size(); ++i) {
        _sumCouplingCells[i]->setMacroscopicMass(0);
        _sumCouplingCells[i]->setMacroscopicMomentum(tarch::la::Vector<dim, double>(0.0));
      }
    }
    for (unsigned int n = 0; n < _totalNumberMDSimulations; ++n) {
      if (0 != _warmupPhase[n])
        continue;
      if (auto* v = dynamic_cast<CouplingCellServiceImpl<LinkedCell, dim>*>(_couplingCellServices[n])) {
        for (unsigned int i = 0; i < v->getCouplingCells().getCouplingCells().size(); ++i) {
          if (_sumCouplingCells.size() <= i) {
            _sumCouplingCells.emplace_back(new coupling::datastructures::CouplingCell<dim>());
          }
          _sumCouplingCells[i]->addMacroscopicMass(v->getCouplingCells().getCouplingCells()[i]->getMacroscopicMass());
          _sumCouplingCells[i]->addMacroscopicMomentum(v->getCouplingCells().getCouplingCells()[i]->getMacroscopicMomentum());
        }
      }
    }
  }

  /** reduce data from MD simulations, averages over them (only macroscopic
   * mass/momentum is considered) and writes the result back into
   * couplingCellsFromMacroscopicSolver. */
  double reduceFromMD2Macro(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver, const I00* const indices) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_NO)
    // Fall back on sequential operation when MPI is not available (avoids redundant implementation)
    return sendFromMD2Macro(couplingCellsFromMacroscopicSolver, indices);
#else
    double res = 0;
    const unsigned int size = (unsigned int)couplingCellsFromMacroscopicSolver.size();

    preprocessingForMD2Macro(indices, size);

    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (nullptr == _couplingCellServices[i] || 0 != _warmupPhase[i])
        continue;
      _couplingCellServices[i]->sendFromMD2MacroPreProcess();
      res += _couplingCellServices[i]->applyFilterPipeline();
    }

    /**
     * local reduction step, has to be executed before EVERY sendFromMD2Macro
     */
    sumUpCouplingCellsFromMamico();

    // reset macroscopic data (only those should be used by macroscopic solver anyway) in cells
    for (coupling::datastructures::CouplingCell<dim>*& couplingCell : _couplingCells) {
      couplingCell->setMacroscopicMass(0.0);
      couplingCell->setMacroscopicMomentum(tarch::la::Vector<dim, double>(0.0));
    }

    std::vector<coupling::sendrecv::DataExchangeFromMD2Macro<dim>*> allDEs(_totalNumberMDSimulations);
    std::vector<std::vector<coupling::datastructures::CouplingCell<dim>*>> allCouplingCellsFromMamico(_totalNumberMDSimulations);
    coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface =
        coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().getMacroscopicSolverInterface();
    unsigned int totalNumberEquilibratedMDSimulations = 0;
    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (nullptr == _couplingCellServices[i] || 0 != _warmupPhase[i])
        continue; // Only reduce using equilibrated MD simulation instances
      allDEs[i] = new coupling::sendrecv::DataExchangeFromMD2Macro<dim>(_macroscopicSolverInterface, _couplingCellServices[i]->getID());
      totalNumberEquilibratedMDSimulations += 1;
    }
    _fromMD2Macro.reduceFromMD2Macro(allDEs, couplingCellsFromMacroscopicSolver, indices, _sumCouplingCells);

    // average data
    for (unsigned int i = 0; i < size; i++) {
      _couplingCells[i]->setMacroscopicMass(couplingCellsFromMacroscopicSolver[i]->getMacroscopicMass() / (double)totalNumberEquilibratedMDSimulations);
      _couplingCells[i]->setMacroscopicMomentum(1.0 / (double)totalNumberEquilibratedMDSimulations *
                                                couplingCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
    }

    for (unsigned int i = 0; i < _totalNumberMDSimulations; i++) {
      if (nullptr == _couplingCellServices[i] || 0 != _warmupPhase[i])
        continue; // Only equilibrated MD simulation instances
      _couplingCellServices[i]->sendFromMD2MacroPostProcess();
    }

    // apply post multi instance FilterPipeline on cell data

#ifdef DEBUG_FILTER_PIPELINE
    std::cout << "FP: Now applying post-multi-instance filter pipeline" << std::endl;
#endif

    (*_postMultiInstanceFilterPipeline)();

    // store data in couplingCellsFromMacroscopicSolver
    for (unsigned int i = 0; i < size; i++) {
      couplingCellsFromMacroscopicSolver[i]->setMacroscopicMass(_couplingCells[i]->getMacroscopicMass());
      couplingCellsFromMacroscopicSolver[i]->setMacroscopicMomentum(_couplingCells[i]->getMacroscopicMomentum());
    }
    return res;
#endif
  }

  /** collects data from MD simulations, averages over them (only macroscopic
   * mass/momentum is considered) and writes the result back into
   * couplingCellsFromMacroscopicSolver. */
  double sendFromMD2Macro(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver, const I00* const indices) {
    double res = 0;
    const unsigned int size = (unsigned int)couplingCellsFromMacroscopicSolver.size();

    preprocessingForMD2Macro(indices, size);

    // reset macroscopic data (only those should be used by macroscopic solver
    // anyway) in duplicate
    for (unsigned int i = 0; i < size; i++) {
      _couplingCells[i]->setMacroscopicMass(0.0);
      _couplingCells[i]->setMacroscopicMomentum(tarch::la::Vector<dim, double>(0.0));
    }

    // receive data from each MD simulation and accumulate information in
    // duplicate
    unsigned int totalNumberEquilibratedMDSimulations = 0;
    for (unsigned int l = 0; l < _totalNumberMDSimulations; l++) {
      // std::cout << "Rank " <<
      //_couplingCellServices[l]->getIndexConversion().getThisRank() << ":
      // Send from MD to Macro for Simulation no. " << l << std::endl;
      if (_couplingCellServices[l] != nullptr && _warmupPhase[l] == 0) {
        res += _couplingCellServices[l]->sendFromMD2Macro(couplingCellsFromMacroscopicSolver, indices);
        for (unsigned int i = 0; i < size; i++) {
          _couplingCells[i]->addMacroscopicMass(couplingCellsFromMacroscopicSolver[i]->getMacroscopicMass());
          _couplingCells[i]->addMacroscopicMomentum(couplingCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
        }
        totalNumberEquilibratedMDSimulations += 1;
      }
    }

    // average data
    for (unsigned int i = 0; i < size; i++) {
      _couplingCells[i]->setMacroscopicMass(_couplingCells[i]->getMacroscopicMass() / (double)totalNumberEquilibratedMDSimulations);
      _couplingCells[i]->setMacroscopicMomentum(1.0 / (double)totalNumberEquilibratedMDSimulations * _couplingCells[i]->getMacroscopicMomentum());
    }

    // apply post multi instance FilterPipeline on cell data

#ifdef DEBUG_FILTER_PIPELINE
    std::cout << "FP: Now applying post-multi-instance filter pipeline" << std::endl;
#endif

    (*_postMultiInstanceFilterPipeline)();

    // store data in couplingCellsFromMacroscopicSolver
    for (unsigned int i = 0; i < size; i++) {
      couplingCellsFromMacroscopicSolver[i]->setMacroscopicMass(_couplingCells[i]->getMacroscopicMass());
      couplingCellsFromMacroscopicSolver[i]->setMacroscopicMomentum(_couplingCells[i]->getMacroscopicMomentum());
    }

    return res;
  }

  /**
   * preprocessing operations for reducing data from MD instances to macro
   */
  void preprocessingForMD2Macro(const I00* const indices, const unsigned int size) {

    /*
     * If this is first coupling step, we must allocate space for the
     * coupling cells we filter and determine averages with.
     */
    if (_couplingCells.empty()) {
      // Allocate & init _couplingCells
      for (unsigned int c = _couplingCells.size(); c < size; c++)
        _couplingCells.push_back(new coupling::datastructures::CouplingCell<dim>());
    }

    /*
     * If this is the first coupling step, copy cells to be in format compatible
     * with filtering system. Can be optimized if indices don't change
     * dynamically.
     */
    if (_cellIndices.empty()) {
      for (unsigned int i = 0; i < size; i++)
        _cellIndices.push_back(indices[i]);
    }

#ifdef ENABLE_POST_MULTI_INSTANCE_FILTERING
    /*
     * If this is the first coupling step, we must init the post multi instance
     * filter pipeline operating on averaged cell data. If you wish to not use post multi-instance filtering
     * in deployment, you can simply leave the corresponding XML-Subtag empty.
     * The ENABLE_POST_MULTI_INSTANCE_FILTERING flag is used for debugging purposes and shall be removed later.
     */
    if (_postMultiInstanceFilterPipeline == nullptr) {
      // Init filter pipeline
      _postMultiInstanceFilterPipeline = new coupling::filtering::FilterPipeline<dim>(_couplingCells, coupling::filtering::Scope::postMultiInstance,
                                                                                      _multiMDService, _filterPipelineConfiguration.c_str());
    }
#endif
  }

  /** removes the last simulation which has been added.
   *
   */
  unsigned int rmMDSimulation(coupling::InstanceHandling<LinkedCell, dim>& instanceHandling, const unsigned int& index) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (_couplingCellServices[index] == nullptr) {
      std::cout << "Rank " << _multiMDService.getGlobalRank() << ": _couplingCellService at " << index << " == NULL" << std::endl;
    }
#endif

    delete _couplingCellServices[index];
    _couplingCellServices[index] = nullptr;

    if (index >= _blockOffset && index < _blockOffset + _localNumberMDSimulations) {
      unsigned int iSim = index - _blockOffset;
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (_multiMDService.getGlobalRank() == 0)
        std::cout << "MultiMDCellService: removing instance " << iSim << "..." << std::endl;
#endif
      instanceHandling.rmMDSimulation(iSim);
    }

    _warmupPhase[index] = -1;

    return index;
  }

  /** In case there is exactly one free slot available per
   *  local comm, remove one slot per comm (that is, one block of simulations)
   *  */
  void removeSimulationBlock() {

    unsigned int newLocalNumberMDSimulations = _localNumberMDSimulations - 1;
    unsigned int newTotalNumberMDSimulations = _totalNumberMDSimulations - _multiMDService.getNumberLocalComms();
    unsigned int newBlockOffset = newLocalNumberMDSimulations * _topologyOffset / _intNumberProcesses;

    auto** newCouplingCellServices = new CouplingCellService<dim>*[newTotalNumberMDSimulations];

    for (unsigned int i = 0; i < _multiMDService.getNumberLocalComms(); ++i) {
      for (unsigned int j = 0; j < newLocalNumberMDSimulations; ++j) {
        unsigned int index = i * _localNumberMDSimulations + j;
        unsigned int newIndex = i * newLocalNumberMDSimulations + j;

        newCouplingCellServices[newIndex] = _couplingCellServices[index];
        if (newIndex < newBlockOffset || newIndex >= newBlockOffset + newLocalNumberMDSimulations) {
          coupling::indexing::IndexingService<dim>::getInstance().updateTopologyOffset((newIndex / newLocalNumberMDSimulations) * _intNumberProcesses);
        }
      }

      auto pos = _warmupPhase.begin() + newLocalNumberMDSimulations * i - i;
      _warmupPhase.erase(pos);
    }

    // Update local variables
    _localNumberMDSimulations = newLocalNumberMDSimulations;
    _totalNumberMDSimulations = newTotalNumberMDSimulations;
    _blockOffset = newBlockOffset;

    delete[] _couplingCellServices;
    _couplingCellServices = newCouplingCellServices;

    //_multiMDService.removeMDSimulationBlock();

    //_listActiveMDSimulations.pop_back();
  }

  /** In order to make space for new simulation slots
   *  We need to add one slot to the local coupling cell service
   * implementations
   *  and also one slot per macroOnly cell service.
   *  As the topologyOffset per MacroOnly cell service changes,
   *  we re-initialize these as well.
   */
  void addSimulationBlock() {

    unsigned int newLocalNumberMDSimulations = _localNumberMDSimulations + 1;
    unsigned int newTotalNumberMDSimulations = _totalNumberMDSimulations + _multiMDService.getNumberLocalComms();
    unsigned int newBlockOffset = newLocalNumberMDSimulations * _topologyOffset / _intNumberProcesses;

    auto** newCouplingCellServices = new CouplingCellService<dim>*[newTotalNumberMDSimulations];
    for (unsigned int i = 0; i < _multiMDService.getNumberLocalComms(); ++i) {
      for (unsigned int j = 0; j < _localNumberMDSimulations; ++j) {
        unsigned int index = i * _localNumberMDSimulations + j;
        unsigned int newIndex = i * newLocalNumberMDSimulations + j;
        newCouplingCellServices[newIndex] = _couplingCellServices[index];
        if (newIndex < newBlockOffset || newIndex >= newBlockOffset + newLocalNumberMDSimulations) {
          coupling::indexing::IndexingService<dim>::getInstance().updateTopologyOffset((newIndex / newLocalNumberMDSimulations) * _intNumberProcesses);
        }
      }
      newCouplingCellServices[(i + 1) * newLocalNumberMDSimulations - 1] = nullptr;

      auto pos = _warmupPhase.begin() + ((i + 1) * newLocalNumberMDSimulations - 1);
      _warmupPhase.insert(pos, 0);
    }

    // Update local variables
    _localNumberMDSimulations = newLocalNumberMDSimulations;
    _totalNumberMDSimulations = newTotalNumberMDSimulations;
    _blockOffset = newBlockOffset;
    delete[] _couplingCellServices;
    _couplingCellServices = newCouplingCellServices;
  }

  /** Adds CouplingCellService at appropriate slot
   *  @return true if this process needs another md simulation initialized
   *          false otherwise
   * */
  unsigned int addMDSimulation(coupling::InstanceHandling<LinkedCell, dim>& instanceHandling,
                               coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface, const unsigned int& slot) {

    if (slot >= _totalNumberMDSimulations) {
      std::cout << "ERROR coupling::services::MultiMDCellService::addMDSimulation(): "
                   "Invalid slot "
                << slot << "!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (_couplingCellServices[slot] != nullptr) {
      std::cout << "ERROR! "
                   "coupling::services::MultiMDCellService::addMDSimulation(): "
                   "Simulation at "
                << slot << " already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    std::stringstream filestem;
    filestem << "restart_checkpoint_" << _multiMDService.getGlobalRank() / computeScalarNumberProcesses();
    instanceHandling.writeCheckpoint(filestem.str().c_str(), 0);

    if (slot < _blockOffset || slot >= _blockOffset + _localNumberMDSimulations) {

      _couplingCellServices[slot] = new coupling::services::CouplingCellServiceDummy<dim>(
          slot, macroscopicSolverInterface, _mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), _multiMDService.getGlobalRank(),
          _mdConfiguration.getDomainConfiguration().getGlobalDomainSize(), _mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
          _mamicoConfiguration.getParallelTopologyConfiguration(), _mamicoConfiguration.getCouplingCellConfiguration(),
          (slot / _localNumberMDSimulations) * _intNumberProcesses);
    } else {

      unsigned localIndex = slot - _blockOffset;
      auto* mdSolverInterface = instanceHandling.addMDSimulation(slot, localIndex);

      _couplingCellServices[slot] = new coupling::services::CouplingCellServiceImpl<LinkedCell, dim>(
          slot, mdSolverInterface, macroscopicSolverInterface, _mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), _multiMDService.getGlobalRank(),
          _mamicoConfiguration.getParticleInsertionConfiguration(), _mamicoConfiguration.getMomentumInsertionConfiguration(),
          _mamicoConfiguration.getBoundaryForceConfiguration(), _mamicoConfiguration.getTransferStrategyConfiguration(),
          _mamicoConfiguration.getParallelTopologyConfiguration(), _mamicoConfiguration.getThermostatConfiguration(),
          _mdConfiguration.getSimulationConfiguration().getNumberOfTimesteps(), _mamicoConfiguration.getCouplingCellConfiguration(),
          _filterPipelineConfiguration.c_str(), _multiMDService, _topologyOffset, _tws);
      instanceHandling.getSimpleMD()[localIndex]->setCouplingCellService((_couplingCellServices[slot]));
      _couplingCellServices[slot]->perturbateVelocity();
      if (_couplingCellServices[slot]->getFilterPipeline() == nullptr) {
        _couplingCellServices[slot]->initFiltering();
      }
    }

    _warmupPhase[slot] = 10;

    return slot;
  }

  unsigned int getLocalNumberOfMDSimulations() const { return _localNumberMDSimulations; }

  void finishCycle(const unsigned int& cycle, coupling::InstanceHandling<LinkedCell, dim>& instanceHandling) {
    // TODO this should be move to somewhere else (MDMediator?)
    for (unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
      if (_warmupPhase[i] > 0) {
        _warmupPhase[i] -= 1;
        if (_warmupPhase[i] == 0 && i >= _blockOffset && i < _blockOffset + _localNumberMDSimulations) {
          instanceHandling.switchOnCoupling(i - _blockOffset);
        }
      }
    }
    // writeCheckpoint(cycle, instanceHandling);
    // TODO call directly instanceHandling.writeCheckpoint();
  }

  void writeCheckpoint(const unsigned int& cycle, const coupling::InstanceHandling<LinkedCell, dim>& instanceHandling) const {
    if (cycle % 10 != 0)
      return;
    std::stringstream filestem;
    filestem << "restart_checkpoint_" << _multiMDService.getGlobalRank() / computeScalarNumberProcesses();
    instanceHandling.writeCheckpoint(filestem.str().c_str(), 0);
  }

  // Must be called after construction of MultiMDCellService, but before the
  // first coupling step.
  void constructFilterPipelines() {
    for (unsigned int md = 0; md < _totalNumberMDSimulations; md++) {
      // get coupling cell service of instance
      auto& mcs = _couplingCellServices[md];

      // only Impl instances of MCS contain filtering
      if (dynamic_cast<coupling::services::CouplingCellServiceImpl<LinkedCell, dim>*>(mcs) == nullptr)
        continue;

      if (mcs->getFilterPipeline() == nullptr) {
        mcs->initFiltering();
      }
    }
  }

private:
  unsigned int computeTopologyOffset() const {
    // determine topology offset of this rank
    const unsigned int intNumberProcesses = computeScalarNumberProcesses();
    const unsigned int topologyOffset = (_multiMDService.getGlobalRank() / intNumberProcesses) * intNumberProcesses;
    return topologyOffset;
  }

  unsigned int computeScalarNumberProcesses() const {
    unsigned int np = _multiMDService.getNumberProcessesPerMDSimulation()[0];
    for (unsigned int d = 1; d < dim; d++) {
      np = np * _multiMDService.getNumberProcessesPerMDSimulation()[d];
    }
    return np;
  }

  coupling::services::CouplingCellService<dim>*
  createCouplingCellServiceDummy(unsigned int ID, coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface,
                                 simplemd::configurations::MolecularDynamicsConfiguration& mdConfiguration, tarch::utils::MultiMDService<dim>& multiMDService,
                                 coupling::configurations::MaMiCoConfiguration<dim>& mamicoConfiguration, unsigned int topologyOffset) const {
    return new coupling::services::CouplingCellServiceDummy<dim>(
        ID, macroscopicSolverInterface, mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), multiMDService.getGlobalRank(),
        mdConfiguration.getDomainConfiguration().getGlobalDomainSize(), mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
        mamicoConfiguration.getParallelTopologyConfiguration(), mamicoConfiguration.getCouplingCellConfiguration(), topologyOffset);
  }

  unsigned int _localNumberMDSimulations;
  /** number of MD simulations run on the current rank. This can differ for
  different blocks, i.e. different topologyOffset values. */
  unsigned int _totalNumberMDSimulations; /** total number of MD simulations */
  coupling::services::CouplingCellService<dim>** _couplingCellServices;
  /** pointers of CouplingCellService type, one for each MD simulation */
  tarch::utils::MultiMDService<dim>& _multiMDService;
  unsigned int _topologyOffset; /** topology offset*/

  const int _tws;
  const unsigned int _intNumberProcesses;

  // const coupling::IndexConversion<dim> _indexConversion; /* Used for index
  // conversions during filtering TODO after merge with dynamic-md*/
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;    /** used to store in MD data in sendFromMD2Macro */
  std::vector<I00> _cellIndices;                                               /** used to store in indexing of the above */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _sumCouplingCells; /** used to reduce all local coupling cells before sending from md 2 macro */

  simplemd::configurations::MolecularDynamicsConfiguration& _mdConfiguration;
  coupling::configurations::MaMiCoConfiguration<dim>& _mamicoConfiguration;
  const std::string _filterPipelineConfiguration;
  coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface;

  unsigned int _blockOffset;
  std::vector<bool> _listActiveMDSimulations; /** One entry per (in-)active md
  simulation, totals to _totalNumberMDSimulations */
  unsigned int _nextFreeBlock;
  /** Points to the next block, to which a simulation should be added. */
  std::vector<int> _warmupPhase;
  /** Counts the number of remaining warmup cycles after this simulation has
               been added.
               This is only relevent for simulations added during the
               simulation.
           */
  /*
   * Analogon to CouplingCellService's FilterPipeline.
   * Is applied during this->sendFromMD2Macro.
   */
  coupling::filtering::FilterPipeline<dim>* _postMultiInstanceFilterPipeline;
  /* Buffer for copying data from MD to macro */
  coupling::sendrecv::FromMD2Macro<coupling::datastructures::CouplingCell<dim>, dim> _fromMD2Macro;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
