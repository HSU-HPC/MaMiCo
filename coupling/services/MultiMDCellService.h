// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_

#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/services/MacroscopicCellService.h"
#include "coupling/services/MacroscopicCellServiceMacroOnly.h"
#include "coupling/interface/MDSimulationFactory.h"

namespace coupling{
  namespace services {
    template<class LinkedCell,unsigned int dim>
    class MultiMDCellService;
  }
}


/** service to couple various MD simulations to a single instance of a macroscopic solver.
 *  We currently consider only one layout: global ranks are linearly split into blocks of size numberProcesses. Each block may run several MD simulations.
 *  Each of these simulations has one mdSolverInterface (handed over in the constructor as vector).
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::services::MultiMDCellService {
  public:

    MultiMDCellService(
      std::vector<coupling::interface::MDSolverInterface<LinkedCell,dim>* > mdSolverInterfaces, // MD solver interfaces for each MD simulation (which uses this rank)
      coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface,
      simplemd::configurations::MolecularDynamicsConfiguration & mdConfiguration,
      coupling::configurations::MaMiCoConfiguration<dim> & mamicoConfiguration,
      const char* filterPipelineConfiguration,
      tarch::utils::MultiMDService<dim> & multiMDService,
      int tws = 0
    ) :
      
      _multiMDService(multiMDService),
      _tws(tws),
      _intNumberProcesses(computeScalarNumberProcesses()),
      _mdConfiguration(mdConfiguration),
      _mamicoConfiguration(mamicoConfiguration),
      _filterPipelineConfiguration(filterPipelineConfiguration),
      _indexConversion(initIndexConversion(
        _mamicoConfiguration.getMacroscopicCellConfiguration().getMacroscopicCellSize(),
        _multiMDService.getNumberProcessesPerMDSimulation(),
        _multiMDService.getRank(),
        _mdConfiguration.getDomainConfiguration().getGlobalDomainSize(),
        _mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
        _mamicoConfiguration.getParallelTopologyConfiguration().getParallelTopologyType(), 
        computeTopologyOffset()
      ))
      {
        _topologyOffset = computeTopologyOffset();
        _localNumberMDSimulations = multiMDService.getLocalNumberOfMDSimulations();
        _totalNumberMDSimulations = multiMDService.getTotalNumberOfMDSimulations();

        _blockOffset = _localNumberMDSimulations*_topologyOffset/_intNumberProcesses;

        _listActiveMDSimulations = std::vector<bool>(_totalNumberMDSimulations, true);
        _nextFreeBlock = _multiMDService.getNumberLocalComms()-1;
        _warmupPhase = std::vector<unsigned int>(_totalNumberMDSimulations, 0);

        const tarch::la::Vector<dim,double> mdDomainSize(mdSolverInterfaces[0]->getGlobalMDDomainSize());
        const tarch::la::Vector<dim,double> mdDomainOffset(mdSolverInterfaces[0]->getGlobalMDDomainOffset());

        // determine globally unique IDs for each macroscopic cell service. Assumptions:
      // - the global parallel topology is subdivided into equally sized blocks of ranks
      // - on each block of ranks, the SAME number of MD simulations is executed. This yields that the variable _topologyOffset together with the local ID of an MD simulation yields a unique global ID
      // - an MD simulation executes on one full block of ranks
      _macroscopicCellServices = new coupling::services::MacroscopicCellService<dim>* [_totalNumberMDSimulations];
      if (_macroscopicCellServices==NULL){std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices==NULL!" << std::endl; exit(EXIT_FAILURE);}

      // allocate all macroscopic cell services for macro-only-solver BEFORE topology offset
      // -> we have localNumberMDSimulations that run per block of ranks on intNumberProcesses
      //    -> this yields localNumberMDSimulations*topologyOffset/intNumberProcesses MD simulations before the actual block of ranks
      for (unsigned int i = 0; i < _blockOffset; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
                                        i, macroscopicSolverInterface, 
                                        mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), 
                                        _multiMDService.getGlobalRank(), 
                                        mdConfiguration.getDomainConfiguration().getGlobalDomainSize(), 
                                        mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
                                        mamicoConfiguration.getParallelTopologyConfiguration(), 
                                        mamicoConfiguration.getMacroscopicCellConfiguration(), (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro- and micro-interactions
      for (unsigned int i = _blockOffset; i<_localNumberMDSimulations+_blockOffset; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>(
                                        i, mdSolverInterfaces[i-_blockOffset], 
                                        macroscopicSolverInterface, 
                                        mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), 
                                        _multiMDService.getGlobalRank(), 
                                        mamicoConfiguration.getParticleInsertionConfiguration(),
                                        mamicoConfiguration.getMomentumInsertionConfiguration(), 
                                        mamicoConfiguration.getBoundaryForceConfiguration(), 
                                        mamicoConfiguration.getTransferStrategyConfiguration(), 
                                        mamicoConfiguration.getNoiseReductionConfiguration(), 
                                        mamicoConfiguration.getParallelTopologyConfiguration(),
                                        mdConfiguration.getSimulationConfiguration().getNumberOfTimesteps(), 
                                        mamicoConfiguration.getMacroscopicCellConfiguration(), 
                                        _filterPipelineConfiguration,
                                        multiMDService, _topologyOffset, _tws
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro-only-solver AFTER topology offset
      for (unsigned int i = _blockOffset+_localNumberMDSimulations; i < _totalNumberMDSimulations; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
                                        i, macroscopicSolverInterface, 
                                        mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), 
                                        _multiMDService.getGlobalRank(), 
                                        mdConfiguration.getDomainConfiguration().getGlobalDomainSize(), 
                                        mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
                                        mamicoConfiguration.getParallelTopologyConfiguration(), 
                                        mamicoConfiguration.getMacroscopicCellConfiguration(), (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }

      _mdConfiguration.getDomainConfigurationNonConst().setInitFromCheckpoint(true);
      std::stringstream filestem;
      filestem << "restart_checkpoint_" << (_multiMDService.getGlobalRank()+1) / _multiMDService.getNumberLocalComms() << "_0";
      _mdConfiguration.getDomainConfigurationNonConst().setCheckpointFilestem(filestem.str());
      _mdConfiguration.getDomainConfigurationNonConst().setInitFromSequentialCheckpoint(false);
      if(_multiMDService.getLocalSize() > 1) {
        //_mdConfiguration.getDomainConfigurationNonConst().setInitFromSequentialCheckpoint(false);
      }
      else {
        //_mdConfiguration.getDomainConfigurationNonConst().setInitFromSequentialCheckpoint(true);
      }
      

    }

    /*MultiMDCellService(
      std::vector<coupling::interface::MDSolverInterface<LinkedCell,dim>* > mdSolverInterfaces, // MD solver interfaces for each MD simulation (which uses this rank)
      coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface,
      tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int rank,
      unsigned int totalNumberMDSimulations,
      const coupling::configurations::ParticleInsertionConfiguration &particleInsertionConfiguration,
      const coupling::configurations::MomentumInsertionConfiguration &momentumInsertionConfiguration,
      const coupling::configurations::BoundaryForceConfiguration<dim> &boundaryForceConfiguration,
      const coupling::configurations::TransferStrategyConfiguration<dim>& transferStrategyConfiguration,
      const coupling::configurations::NoiseReductionConfiguration &noiseReductionConfiguration,
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,
      unsigned int numberMDTimestepsPerCouplingCycle,
      const coupling::configurations::MacroscopicCellConfiguration<dim> &macroscopicCellConfiguration,
	    const char* filterPipelineConfiguration,
       tarch::utils::MultiMDService<dim>& multiMDService, int tws = 0
    ): _localNumberMDSimulations((unsigned int)mdSolverInterfaces.size()), _totalNumberMDSimulations(totalNumberMDSimulations),
       _macroscopicCellServices(NULL), 
       _multiMDService(multiMDService), 
       _tws(tws), 
       _intNumberProcesses(computeScalarNumberProcesses())
       {

      _topologyOffset = computeTopologyOffset();

      const tarch::la::Vector<dim,double> mdDomainSize(mdSolverInterfaces[0]->getGlobalMDDomainSize());
      const tarch::la::Vector<dim,double> mdDomainOffset(mdSolverInterfaces[0]->getGlobalMDDomainOffset());

      // determine globally unique IDs for each macroscopic cell service. Assumptions:
      // - the global parallel topology is subdivided into equally sized blocks of ranks
      // - on each block of ranks, the SAME number of MD simulations is executed. This yields that the variable _topologyOffset together with the local ID of an MD simulation yields a unique global ID
      // - an MD simulation executes on one full block of ranks
      _macroscopicCellServices = new coupling::services::MacroscopicCellService<dim>* [_totalNumberMDSimulations];
      if (_macroscopicCellServices==NULL){std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices==NULL!" << std::endl; exit(EXIT_FAILURE);}

      // allocate all macroscopic cell services for macro-only-solver BEFORE topology offset
      // -> we have localNumberMDSimulations that run per block of ranks on intNumberProcesses
      //    -> this yields localNumberMDSimulations*topologyOffset/intNumberProcesses MD simulations before the actual block of ranks
      for (unsigned int i = 0; i < _localNumberMDSimulations*_topologyOffset/_intNumberProcesses; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
                                        i, macroscopicSolverInterface, numberProcesses, _multiMDService.getGlobalRank(), mdDomainSize, mdDomainOffset,
                                        parallelTopologyConfiguration, macroscopicCellConfiguration, (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro- and micro-interactions
      for (unsigned int i = _localNumberMDSimulations*_topologyOffset/_intNumberProcesses; i<_localNumberMDSimulations*(_topologyOffset+_intNumberProcesses)/_intNumberProcesses; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>(
                                        i, mdSolverInterfaces[i-_localNumberMDSimulations*_topologyOffset/_intNumberProcesses], macroscopicSolverInterface, numberProcesses, _multiMDService.getGlobalRank(), particleInsertionConfiguration,
                                        momentumInsertionConfiguration, boundaryForceConfiguration, transferStrategyConfiguration, noiseReductionConfiguration, parallelTopologyConfiguration,
                                        numberMDTimestepsPerCouplingCycle, macroscopicCellConfiguration, filterPipelineConfiguration, multiMDService, _topologyOffset, _tws
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro-only-solver AFTER topology offset
      for (unsigned int i = _localNumberMDSimulations*(_topologyOffset+_intNumberProcesses)/_intNumberProcesses; i < _totalNumberMDSimulations; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
                                        i, macroscopicSolverInterface, numberProcesses, _multiMDService.getGlobalRank(), mdDomainSize, mdDomainOffset,
                                        parallelTopologyConfiguration, macroscopicCellConfiguration, (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
    }*/

    ~MultiMDCellService(){
      for (unsigned int i = 0; i < _totalNumberMDSimulations; i++){
        if (_macroscopicCellServices[i]!=nullptr){
          delete _macroscopicCellServices[i]; 
          _macroscopicCellServices[i] = nullptr; 
          }
      }
      if (_macroscopicCellServices!=NULL){ 
        delete [] _macroscopicCellServices; 
        _macroscopicCellServices=NULL; 
      }
      if (_indexConversion!=nullptr) {
        delete _indexConversion;
        _indexConversion = nullptr;
      }
    }

    /** get access to macroscopic cell services. This is required potentially by each MD simulation to incorporate cell-local thermostat, mass and momentum transfer */
    coupling::services::MacroscopicCellService<dim>& getMacroscopicCellService(unsigned int localIndex){
      if (localIndex<_localNumberMDSimulations){ return *(_macroscopicCellServices[_topologyOffset*_localNumberMDSimulations/_intNumberProcesses+localIndex]); }
      else{ std::cout << "ERROR MultiMDCellService::getMacroscopicCellService(localIndex): localIndex >_localNumberMDSimulations-1!" << std::endl; exit(EXIT_FAILURE);}
    }

    /** Forward Plotting to macroscopic cell services if not NULL! */
    void plotEveryMacroscopicTimestepforMacroscopicCellService(unsigned int localIndex, int cycle) {
      if (_macroscopicCellServices[_topologyOffset*_localNumberMDSimulations/_intNumberProcesses+localIndex] == nullptr) return;
      _macroscopicCellServices[_topologyOffset*_localNumberMDSimulations/_intNumberProcesses+localIndex]->plotEveryMacroscopicTimestep(cycle);
    }

    void plotEveryMacroscopicTimestep(int cycle) {
      for(unsigned int i = _blockOffset;i<_blockOffset+_localNumberMDSimulations;++i) {
        if (_macroscopicCellServices[i] == nullptr) continue;
        _macroscopicCellServices[i]->plotEveryMacroscopicTimestep(cycle);
      }
    }

    void computeAndStoreTemperature(double temp) {
      for(unsigned int i = _blockOffset;i<_blockOffset+_localNumberMDSimulations;++i) {
        if (_macroscopicCellServices[i] == nullptr) continue;
        _macroscopicCellServices[i]->computeAndStoreTemperature(temp);
      }
    }


    /** broadcasts data from macroscopic solver to all MD simulations */
    void sendFromMacro2MD(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ){
      // just send information to all MD instances. This is currently sequentialized inside MacroscopicCellService/SendRecvBuffer
      for (unsigned int i = 0; i < _totalNumberMDSimulations; i++){
        //std::cout << "Rank " << _macroscopicCellServices[i]->getIndexConversion().getThisRank() << ": Send from macro to MD for Simulation no. " << i << std::endl;
        if(_macroscopicCellServices[i] != nullptr) {
          _macroscopicCellServices[i]->sendFromMacro2MD(macroscopicCellsFromMacroscopicSolver,globalCellIndicesFromMacroscopicSolver);
        }
      }
    }

    /** collects data from MD simulations, averages over them (only macroscopic mass/momentum is considered) and writes the result back into macroscopicCellsFromMacroscopicSolver. */
    double sendFromMD2Macro(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ){
      double res = 0;
      // duplicate macroscopicCells
      const unsigned int size = (unsigned int) macroscopicCellsFromMacroscopicSolver.size();
      coupling::datastructures::MacroscopicCell<dim>* duplicate = new coupling::datastructures::MacroscopicCell<dim>[size];
      if (duplicate==NULL){std::cout << "ERROR coupling::services::MultiMDCellService::sendFromMD2Macro(...): duplicate==NULL!" << std::endl; exit(EXIT_FAILURE);}

      // reset macroscopic data (only those should be used by macroscopic solver anyway) in duplicate
      for (unsigned int i = 0; i < size; i++){
        duplicate[i].setMacroscopicMass(0.0); duplicate[i].setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
      }

      // receive data from each MD simulation and accumulate information in duplicate
      for (unsigned int l = 0; l < _totalNumberMDSimulations; l++){
          //std::cout << "Rank " << _macroscopicCellServices[l]->getIndexConversion().getThisRank() << ": Send from MD to Macro for Simulation no. " << l << std::endl;
          if (_macroscopicCellServices[l] != nullptr &&  _warmupPhase[l] > 0)  {
          res += _macroscopicCellServices[l]->sendFromMD2Macro(macroscopicCellsFromMacroscopicSolver,globalCellIndicesFromMacroscopicSolver);
          for (unsigned int i = 0; i < size; i++){
            duplicate[i].addMacroscopicMass(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMass());
            duplicate[i].addMacroscopicMomentum(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
          }
        }
      }
      // average data and store it in macroscopicCellsFromMacroscopicSolver
      for (unsigned int i = 0; i < size; i++){
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMass(duplicate[i].getMacroscopicMass()/_totalNumberMDSimulations);
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMomentum((1.0/_totalNumberMDSimulations)*duplicate[i].getMacroscopicMomentum());
      }

      // free duplicate
      delete [] duplicate;
      return res;
    }
/*
    auto& getIndexConversion() { 
      for(unsigned int i=0;i<_totalNumberMDSimulations;++i) {
        if (_macroscopicCellServices[i] != NULL) {
          return _macroscopicCellServices[i]->getIndexConversion();
        }
      }
      std::cout << "ERROR coupling::services::MultiMDCellService::getIndexConversion: No valid MacroscopicCellService found!" << std::endl;
      exit(EXIT_FAILURE);
    }
    */


    /** removes the last simulation which has been added.
     *  
     */
    unsigned int rmMDSimulation(std::vector<coupling::interface::MDSolverInterface<LinkedCell,dim>* > & mdSolverInterfaces
                         , std::vector<coupling::interface::MDSimulation*> & simpleMD) {
      if(_localNumberMDSimulations < 2) {
        std::cout << "INFO MultiMDCellService::rmMDSimulation() : Cannot remove MD simulation, only one is left!" << std::endl;
        return 0;
      }
      unsigned int index = getLastReservedSlot();

      delete _macroscopicCellServices[index];
      _macroscopicCellServices[index] = nullptr;

      
      if(index >= _blockOffset && index < _blockOffset + _localNumberMDSimulations) {
        unsigned int iSim = index - _blockOffset;
        simpleMD[iSim]->shutdown();
        delete simpleMD[iSim];
        simpleMD[iSim] = nullptr;
        simpleMD.erase(simpleMD.begin()+iSim);
        delete mdSolverInterfaces[iSim];
        mdSolverInterfaces[iSim] = nullptr;
        mdSolverInterfaces.erase(mdSolverInterfaces.begin()+iSim);
      }

      /** Check if there is a free block of simulations,
       *  and, if so, remove it. Then reset _nextFreeBlock.
       * */
      if(index == _totalNumberMDSimulations-1) {
        removeSimulationBlock();
        _nextFreeBlock = _multiMDService.getNumberLocalComms()-1;
      } 

      return index;
    }

    /** In case there is exactly one free slot available per
     *  local comm, remove one slot per comm (that is, one block of simulations)
     *  */
    void removeSimulationBlock() {

      unsigned int newLocalNumberMDSimulations = _localNumberMDSimulations - 1;
      unsigned int newTotalNumberMDSimulations = _totalNumberMDSimulations - _multiMDService.getNumberLocalComms();
      unsigned int newBlockOffset = newLocalNumberMDSimulations*_topologyOffset/_intNumberProcesses;

      auto ** newMacroscopicCellServices = new MacroscopicCellService<dim> * [newTotalNumberMDSimulations];

      for(unsigned int i=0;i<_multiMDService.getNumberLocalComms();++i) {
        for(unsigned int j=0;j<newLocalNumberMDSimulations;++j) {
          unsigned int index = i * _localNumberMDSimulations + j;
          unsigned int newIndex = i * newLocalNumberMDSimulations + j;

          newMacroscopicCellServices[newIndex] = _macroscopicCellServices[index];
          if (newIndex < newBlockOffset || newIndex >= newBlockOffset + newLocalNumberMDSimulations) {
            // Need to update topology offset in indexconversion->paralleltopology
            newMacroscopicCellServices[newIndex]->updateIndexConversion((newIndex/newLocalNumberMDSimulations) * _intNumberProcesses);
          }
        }

        auto pos = _warmupPhase.begin() + newLocalNumberMDSimulations * i - i;
        _warmupPhase.erase(pos);
      }

      // Update local variables
      _localNumberMDSimulations = newLocalNumberMDSimulations;
      _totalNumberMDSimulations = newTotalNumberMDSimulations;
      _blockOffset = newBlockOffset;

      delete [] _macroscopicCellServices;
      _macroscopicCellServices = newMacroscopicCellServices;

      _multiMDService.removeMDSimulationBlock();

      _listActiveMDSimulations.pop_back();
    }

    /** In order to make space for new simulation slots
     *  We need to add one slot to the local macroscopic cell service implementations
     *  and also one slot per macroOnly cell service.
     *  As the topologyOffset per MacroOnly cell service changes,
     *  we re-initialize these as well.
     */
    void addSimulationBlock() {

      unsigned int newLocalNumberMDSimulations = _localNumberMDSimulations + 1;
      unsigned int newTotalNumberMDSimulations = _totalNumberMDSimulations + _multiMDService.getNumberLocalComms();
      unsigned int newBlockOffset = newLocalNumberMDSimulations*_topologyOffset/_intNumberProcesses;

      auto ** newMacroscopicCellServices = new MacroscopicCellService<dim> * [newTotalNumberMDSimulations];
      for(unsigned int i=0;i<_multiMDService.getNumberLocalComms();++i) {
        for(unsigned int j=0;j<_localNumberMDSimulations;++j) {
          unsigned int index = i * _localNumberMDSimulations + j;
          unsigned int newIndex = i * newLocalNumberMDSimulations + j;
          newMacroscopicCellServices[newIndex] = _macroscopicCellServices[index];
          if (newIndex < newBlockOffset || newIndex >= newBlockOffset + newLocalNumberMDSimulations) {
            // Need to update topologoy offset in indexconversion->paralleltopology
            newMacroscopicCellServices[newIndex]->updateIndexConversion((newIndex/newLocalNumberMDSimulations) * _intNumberProcesses);
          }
        }
        newMacroscopicCellServices[(i+1)*newLocalNumberMDSimulations-1] = nullptr;

        auto pos = _warmupPhase.begin() + ((i+1) * newLocalNumberMDSimulations - 1);
        _warmupPhase.insert(pos, 0);
      }

      // Update local variables
      _localNumberMDSimulations = newLocalNumberMDSimulations;
      _totalNumberMDSimulations = newTotalNumberMDSimulations;
      _blockOffset = newBlockOffset;
      delete [] _macroscopicCellServices;
      _macroscopicCellServices = newMacroscopicCellServices;

      _multiMDService.addMDSimulationBlock();

      _listActiveMDSimulations.push_back(false);
    }

    /** Find and reserve the next available slot
     *  in a round robin fashion.
     */
    unsigned int reserveNextSlot() {
      if(_nextFreeBlock == _multiMDService.getNumberLocalComms()-1) {
        /** If the last communicator gets a new simulation slot,
         *  we first have to add a new block of simulations.
         *  That is, each communicator is appended another slot to hold
         *  a simulation.
         * */
        addSimulationBlock();
      }
      
      unsigned int thisFreeBlock = _nextFreeBlock;
      if(_nextFreeBlock == 0) {
        /** If a simulation slot is reserved on communicator 0,
         *  then we start again at the last one.
         *  The next time a slot is being reserved, we first
         *  have to add a new block of free slots (see condition above)
         *  */
        _nextFreeBlock = _multiMDService.getNumberLocalComms()-1;
      } else {
        _nextFreeBlock -= 1;
      }
      return thisFreeBlock * _localNumberMDSimulations + _localNumberMDSimulations -1;
    }


    /** Find, which slot has been reserved last */
    unsigned int getLastReservedSlot() {
      if(_nextFreeBlock == _multiMDService.getNumberLocalComms()-1) {
        /** In this case, there no actually free slots available
         *  (we first had to add another block)
         *  We thus want to free a slot on the first communicator.
         */
        _nextFreeBlock = 0;
        return _localNumberMDSimulations - 1;       
      }
      unsigned int _thisFreeBlock = _nextFreeBlock;
      _nextFreeBlock += 1;
      return (_thisFreeBlock+1) * _localNumberMDSimulations + _localNumberMDSimulations - 1;
    }


     /** Adds MacroscopicCellService at appropriate slot
      *  @return true if this process needs another md simulation initialized
      *          false otherwise
      * */
     unsigned int addMDSimulation(coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface,
                         std::vector<coupling::interface::MDSolverInterface<LinkedCell,dim>* > & mdSolverInterfaces
                         , std::vector<coupling::interface::MDSimulation*> & simpleMD) {
      
      //TODO find blocks with lowest number of activs Simulations
      //      then add sim to block of highest ID.
      unsigned int slot = reserveNextSlot();

      if(_macroscopicCellServices[slot] != nullptr) {
        std::cout << "ERROR! coupling::services::MultiMDCellService::addMDSimulation(): Simulation at " << slot << " already exists!" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      _listActiveMDSimulations[slot] = true;

      if(slot < _blockOffset || slot >= _blockOffset + _localNumberMDSimulations) {

        _macroscopicCellServices[slot] = 
          new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
            slot, macroscopicSolverInterface, 
            _mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), 
            _multiMDService.getGlobalRank(), 
            _mdConfiguration.getDomainConfiguration().getGlobalDomainSize(), 
            _mdConfiguration.getDomainConfiguration().getGlobalDomainOffset(),
            _mamicoConfiguration.getParallelTopologyConfiguration(), 
            _mamicoConfiguration.getMacroscopicCellConfiguration(), 
            (slot/_localNumberMDSimulations)*_intNumberProcesses
          );
      }
      else {
        auto * mdSim = coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
          _mdConfiguration,_mamicoConfiguration
          #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
          , _multiMDService.getLocalCommunicator()
          #endif
        );

        if(mdSim == NULL) {
          std::cout << "ERROR! coupling::tests::DynamicMDTest::addMDSimulation(): mdSim == NULL!" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        mdSim->init(_multiMDService, slot ); // TODO update variables in multimdservice!!!

        simpleMD.push_back(mdSim);

        mdSolverInterfaces.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().
                                      getMDSolverInterface(_mdConfiguration, _mamicoConfiguration, simpleMD[simpleMD.size()-1]));

        _macroscopicCellServices[slot] = 
          new coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>(
            slot, mdSolverInterfaces[slot-_blockOffset], 
            macroscopicSolverInterface, 
            _mdConfiguration.getMPIConfiguration().getNumberOfProcesses(), 
            _multiMDService.getGlobalRank(), 
            _mamicoConfiguration.getParticleInsertionConfiguration(),
            _mamicoConfiguration.getMomentumInsertionConfiguration(), 
            _mamicoConfiguration.getBoundaryForceConfiguration(), 
            _mamicoConfiguration.getTransferStrategyConfiguration(), 
            _mamicoConfiguration.getNoiseReductionConfiguration(), 
            _mamicoConfiguration.getParallelTopologyConfiguration(),
            _mdConfiguration.getSimulationConfiguration().getNumberOfTimesteps(), 
            _mamicoConfiguration.getMacroscopicCellConfiguration(),
            _filterPipelineConfiguration, 
            _multiMDService, _topologyOffset, _tws
          );
        simpleMD[simpleMD.size()-1]->setMacroscopicCellService((_macroscopicCellServices[slot]));
      }

      _warmupPhase[slot] = 2;

      return slot;
    }

    unsigned int getLocalNumberOfMDSimulations() const { return _localNumberMDSimulations; }

    coupling::IndexConversion<dim> & getIndexConversion() const { return *_indexConversion; }

    void finishCycle(const unsigned int & cycle, const std::vector<coupling::interface::MDSimulation*> simpleMD) {
      for(auto& phaseI : _warmupPhase) {
        if(phaseI > 0) phaseI -= 1;
      }
      writeCheckpoint(cycle, simpleMD);
    }

    void writeCheckpoint(const unsigned int & cycle, 
                          const std::vector<coupling::interface::MDSimulation*> simpleMD) {
      std::stringstream filestem;
      filestem << "restart_checkpoint_" << _multiMDService.getGlobalRank() / _multiMDService.getNumberLocalComms();
      simpleMD[0]->writeCheckpoint(filestem.str().c_str(), 0);
    }
    
  private:
    unsigned int computeTopologyOffset() const {
      // determine topology offset of this rank
      const unsigned int intNumberProcesses = computeScalarNumberProcesses();
      const unsigned int topologyOffset = (_multiMDService.getGlobalRank()/intNumberProcesses)*intNumberProcesses;
      return topologyOffset;
    }

    unsigned int computeScalarNumberProcesses() const {
      unsigned int np = _multiMDService.getNumberProcessesPerMDSimulation()[0];
      for (unsigned int d = 1; d < dim; d++){ np = np*_multiMDService.getNumberProcessesPerMDSimulation()[d]; }
      return np;
    }
    
    coupling::IndexConversion<dim> *
    initIndexConversion(
      tarch::la::Vector<dim,double> macroscopicCellSize,
      tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int rank,
      tarch::la::Vector<dim,double> globalMDDomainSize,
      tarch::la::Vector<dim,double> globalMDDomainOffset,
      coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
      unsigned int topologyOffset) const {

        tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells(0);
        for (unsigned int d = 0; d < dim; d++){
          globalNumberMacroscopicCells[d] = (unsigned int) floor( globalMDDomainSize[d]/macroscopicCellSize[d] + 0.5 );
          if ( fabs(globalNumberMacroscopicCells[d]*macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13 ){
            std::cout << "coupling::services::MultiMDCellService::initIndexConversion(): Deviation of domain size > 1e-13!" << std::endl;
          }
        }
        coupling::IndexConversion<dim> *ic = new coupling::IndexConversion<dim>(globalNumberMacroscopicCells,numberProcesses,rank,globalMDDomainSize,globalMDDomainOffset,parallelTopologyType,topologyOffset);
        if (ic==NULL){std::cout << "coupling::services::MultiMDCellService::initIndexConversion(): ic==NULL!" << std::endl; exit(EXIT_FAILURE);}
        return ic;
    }
    

    unsigned int _localNumberMDSimulations; /** number of MD simulations run on the current rank. This can differ for different blocks, i.e. different topologyOffset values. */
    unsigned int _totalNumberMDSimulations; /** total number of MD simulations */
    coupling::services::MacroscopicCellService<dim> **_macroscopicCellServices; /** pointers of MacroscopicCellService type, one for each MD simulation */
    tarch::utils::MultiMDService<dim>& _multiMDService;
    unsigned int _topologyOffset; /** topology offset*/
    const int _tws;
    const unsigned int _intNumberProcesses;
    
    simplemd::configurations::MolecularDynamicsConfiguration & _mdConfiguration;
    coupling::configurations::MaMiCoConfiguration<dim> & _mamicoConfiguration;
    const char* _filterPipelineConfiguration;

    coupling::IndexConversion<dim> * _indexConversion;

    unsigned int _blockOffset;
    std::vector<bool> _listActiveMDSimulations; /** One entry per (in-)active md simulation, totals to _totalNumberMDSimulations */
    unsigned int _nextFreeBlock; /** Points to the next block, to which a simulation should be added. */ 
    std::vector<unsigned int> _warmupPhase; /** Counts the number of remaining warmup cycles after this simulation has been added.
                                                This is only relevent for simulations added during the simulation.
                                            */

};
#endif // _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
