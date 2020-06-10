// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_

#include "coupling/services/MacroscopicCellService.h"
#include "coupling/services/MacroscopicCellServiceMacroOnly.h"

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
      const tarch::utils::MultiMDService<dim>& multiMDService, int tws = 0
    ): _localNumberMDSimulations((unsigned int)mdSolverInterfaces.size()), _totalNumberMDSimulations(totalNumberMDSimulations),
       _macroscopicCellServices(NULL), _topologyOffset(computeTopologyOffset(numberProcesses,rank)), _tws(tws), _intNumberProcesses(computeScalarNumberProcesses(numberProcesses)) {

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
                                        i, macroscopicSolverInterface, numberProcesses, rank, mdDomainSize, mdDomainOffset,
                                        parallelTopologyConfiguration, macroscopicCellConfiguration, (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro- and micro-interactions
      for (unsigned int i = _localNumberMDSimulations*_topologyOffset/_intNumberProcesses; i<_localNumberMDSimulations*(_topologyOffset+_intNumberProcesses)/_intNumberProcesses; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>(
                                        i, mdSolverInterfaces[i-_localNumberMDSimulations*_topologyOffset/_intNumberProcesses], macroscopicSolverInterface, numberProcesses, rank, particleInsertionConfiguration,
                                        momentumInsertionConfiguration, boundaryForceConfiguration, transferStrategyConfiguration, noiseReductionConfiguration, parallelTopologyConfiguration,
                                        numberMDTimestepsPerCouplingCycle, macroscopicCellConfiguration, multiMDService, _topologyOffset, _tws
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
      // allocate all macroscopic cell services for macro-only-solver AFTER topology offset
      for (unsigned int i = _localNumberMDSimulations*(_topologyOffset+_intNumberProcesses)/_intNumberProcesses; i < _totalNumberMDSimulations; i++){
        _macroscopicCellServices[i] = new coupling::services::MacroscopicCellServiceMacroOnly<dim>(
                                        i, macroscopicSolverInterface, numberProcesses, rank, mdDomainSize, mdDomainOffset,
                                        parallelTopologyConfiguration, macroscopicCellConfiguration, (i/_localNumberMDSimulations)*_intNumberProcesses
                                      );
        if (_macroscopicCellServices[i]==NULL){
          std::cout << "ERROR coupling::services::MultiMDCellService::MultiMDCellService(...): _macroscopicCellServices[" << i << "]==NULL!" << std::endl; exit(EXIT_FAILURE);
        }
      }
    }

    ~MultiMDCellService(){
      for (unsigned int i = 0; i < _totalNumberMDSimulations; i++){
        if (_macroscopicCellServices[i]!=NULL){delete _macroscopicCellServices[i]; _macroscopicCellServices[i] = NULL; }
      }
      if (_macroscopicCellServices!=NULL){ delete [] _macroscopicCellServices; _macroscopicCellServices=NULL; }
    }

    /** get access to macroscopic cell services. This is required potentially by each MD simulation to incorporate cell-local thermostat, mass and momentum transfer */
    coupling::services::MacroscopicCellService<dim>& getMacroscopicCellService(unsigned int localIndex){
      if (localIndex<_localNumberMDSimulations){ return *(_macroscopicCellServices[_topologyOffset*_localNumberMDSimulations/_intNumberProcesses+localIndex]); }
      else{ std::cout << "ERROR MultiMDCellService::getMacroscopicCellService(localIndex): localIndex >_localNumberMDSimulations-1!" << std::endl; exit(EXIT_FAILURE);}
    }

    /** broadcasts data from macroscopic solver to all MD simulations */
    void sendFromMacro2MD(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ){
      // just send information to all MD instances. This is currently sequentialized inside MacroscopicCellService/SendRecvBuffer
      for (unsigned int i = 0; i < _totalNumberMDSimulations; i++){
        //std::cout << "Rank " << _macroscopicCellServices[i]->getIndexConversion().getThisRank() << ": Send from macro to MD for Simulation no. " << i << std::endl;
        _macroscopicCellServices[i]->sendFromMacro2MD(macroscopicCellsFromMacroscopicSolver,globalCellIndicesFromMacroscopicSolver);
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
        res += _macroscopicCellServices[l]->sendFromMD2Macro(macroscopicCellsFromMacroscopicSolver,globalCellIndicesFromMacroscopicSolver);
        for (unsigned int i = 0; i < size; i++){
          duplicate[i].addMacroscopicMass(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMass());
          duplicate[i].addMacroscopicMomentum(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
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

    /** removes select MD simulation */
    void rmMDSimulation(const unsigned int & localIndex) {
      coupling::services::MacroscopicCellService<dim> **mCSTemp = new coupling::services::MacroscopicCellService<dim>* [_totalNumberMDSimulations-1];

      unsigned int i;
      for(i = 0;i<localIndex; ++i) {
        mCSTemp[i] = _macroscopicCellServices[i];
      }
      for(i = localIndex+1;i<_totalNumberMDSimulations;++i) {
        mCSTemp[i] = _macroscopicCellServices[i];
      }

      if(_macroscopicCellServices[localIndex]!=NULL) { 
        delete _macroscopicCellServices[localIndex];
        _macroscopicCellServices[localIndex] = NULL;
      }

      _macroscopicCellServices = mCSTemp;

      _localNumberMDSimulations -= 1;
      _totalNumberMDSimulations -= 1;
    }
    
  private:
    unsigned int computeTopologyOffset(tarch::la::Vector<dim,unsigned int> numberProcesses,unsigned int rank) const {
      // determine topology offset of this rank
      const unsigned int intNumberProcesses = computeScalarNumberProcesses(numberProcesses);
      const unsigned int topologyOffset = (rank/intNumberProcesses)*intNumberProcesses;
      return topologyOffset;
    }

    unsigned int computeScalarNumberProcesses(tarch::la::Vector<dim,unsigned int> numberProcesses) const {
      unsigned int np = numberProcesses[0];
      for (unsigned int d = 1; d < dim; d++){ np = np*numberProcesses[d]; }
      return np;
    }

    unsigned int _localNumberMDSimulations; /** number of MD simulations run on the current rank. This can differ for different blocks, i.e. different topologyOffset values. */
    unsigned int _totalNumberMDSimulations; /** total number of MD simulations */
    coupling::services::MacroscopicCellService<dim> **_macroscopicCellServices; /** pointers of MacroscopicCellService type, one for each MD simulation */
    const unsigned int _topologyOffset; /** topology offset*/
    const int _tws;
    const unsigned int _intNumberProcesses;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
