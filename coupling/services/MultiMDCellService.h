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
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,
      unsigned int numberMDTimestepsPerCouplingCycle,
      const coupling::configurations::MacroscopicCellConfiguration<dim> &macroscopicCellConfiguration,
      const char* filterPipelineConfiguration, //location of the .xml file containing the static filter pipeline config
      const tarch::utils::MultiMDService<dim>& multiMDService, int tws = 0
    ): 
       _localNumberMDSimulations((unsigned int)mdSolverInterfaces.size()), _totalNumberMDSimulations(totalNumberMDSimulations),
       _macroscopicCellServices(NULL), _topologyOffset(computeTopologyOffset(numberProcesses,rank)), _tws(tws), _intNumberProcesses(computeScalarNumberProcesses(numberProcesses)),
       _macroscopicSolverInterface(macroscopicSolverInterface),
       _filterPipelineConfiguration(filterPipelineConfiguration),
       _multiMDService(multiMDService),
       _rank(rank),
       _postMultiInstanceFilterPipeline(nullptr) {

      //If we allow for zero MD instances on ranks, this initializion process would segfault...
      if(mdSolverInterfaces.size() == 0) {
        std::cout << "ERROR: Zero MD instances on rank " << rank << ". Maybe you forgot to adjust number-md-simulations in couette.xml?" << std::endl;
        exit(EXIT_FAILURE);
      }

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
                                        momentumInsertionConfiguration, boundaryForceConfiguration, transferStrategyConfiguration, parallelTopologyConfiguration,
                                        numberMDTimestepsPerCouplingCycle, macroscopicCellConfiguration, filterPipelineConfiguration, multiMDService, _topologyOffset, _tws
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

      //TODO: fix free bug/possible memory leak here
      //delete [] _macroscopicCells.data();
      if (_postMultiInstanceFilterPipeline!=nullptr) delete _postMultiInstanceFilterPipeline;
    }

    /** get access to macroscopic cell services. This is required potentially by each MD simulation to incorporate cell-local thermostat, mass and momentum transfer */
    coupling::services::MacroscopicCellService<dim>& getMacroscopicCellService(unsigned int localIndex){
      if (localIndex<_localNumberMDSimulations){ return *(_macroscopicCellServices[_topologyOffset*_localNumberMDSimulations/_intNumberProcesses+localIndex]); }
      else{ std::cout << "ERROR MultiMDCellService::getMacroscopicCellService(localIndex): localIndex >_localNumberMDSimulations-1!" << std::endl; exit(EXIT_FAILURE);}
    }

    coupling::FilterPipeline<dim>& getPostMultiInstanceFilterPipeline() { return _postMultiInstanceFilterPipeline; }

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

    void bcastFromMacro2MD(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ){

      std::vector<coupling::sendrecv::DataExchangeFromMacro2MD<dim> * > allDEs(_totalNumberMDSimulations);
      std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim> *> > allMacroscopicCellsFromMamico(_totalNumberMDSimulations);
      for(unsigned int i=0;i<_totalNumberMDSimulations;++i) {
        allDEs[i] = new coupling::sendrecv::DataExchangeFromMacro2MD<dim> (
          _macroscopicSolverInterface, _macroscopicCellServices[i]->getIndexConversion(),
          _macroscopicCellServices[i]->getID()
        );
        if(auto * v = dynamic_cast<MacroscopicCellServiceImpl<LinkedCell, dim>*>(_macroscopicCellServices[i])) {
          allMacroscopicCellsFromMamico[i] = v->getMacroscopicCells().getMacroscopicCells();
        }
      }

      for(unsigned int i=0;i<_totalNumberMDSimulations;++i) {
        _macroscopicCellServices[i]->sendFromMacro2MDPreProcess();
      }

      coupling::sendrecv::FromMacro2MD<coupling::datastructures::MacroscopicCell<dim>, dim > fromMacro2MD;
      fromMacro2MD.bcastFromMacro2MD(
        allDEs, macroscopicCellsFromMacroscopicSolver, globalCellIndicesFromMacroscopicSolver, allMacroscopicCellsFromMamico
      );

      for(unsigned int i=0;i<_totalNumberMDSimulations;++i) {
        _macroscopicCellServices[i]->sendFromMacro2MDPostProcess();
      }
      for(coupling::sendrecv::DataExchangeFromMacro2MD<dim> * & de : allDEs) {
        delete de;
        de = nullptr;
      }
    }

    /** Creates the sum over all instances' macroscopic cells, in order to reduce the amount of communication needed.
     */
    void sumUpMacroscopicCellsFromMamico() {
      for(unsigned int n=0;n<_totalNumberMDSimulations;++n) {
        if(auto * v = dynamic_cast<MacroscopicCellServiceImpl<LinkedCell, dim>*>(_macroscopicCellServices[n])) {
          for(unsigned int i = 0; i < v->getMacroscopicCells().getMacroscopicCells().size(); ++i) {
            if(_sumMacroscopicCells.size() <= i) {
              _sumMacroscopicCells.emplace_back(new coupling::datastructures::MacroscopicCell<dim>());
            }
              _sumMacroscopicCells[i]->addMacroscopicMass(v->getMacroscopicCells().getMacroscopicCells()[i]->getMacroscopicMass());
              _sumMacroscopicCells[i]->addMacroscopicMomentum(v->getMacroscopicCells().getMacroscopicCells()[i]->getMacroscopicMomentum());
          }
        }
      }
    }

    double reduceFromMD2Macro(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ) {
      auto * debugTMP = &dynamic_cast<coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>* >(_macroscopicCellServices[0])->getMacroscopicCells();
      unsigned int size = macroscopicCellsFromMacroscopicSolver.size(); // we assume globalCellIndicesFromMacroscopicSolver to be of identical size

      for(unsigned int i = 0; i < _totalNumberMDSimulations; ++i) {
        _macroscopicCellServices[i]->sendFromMD2MacroPreProcess();
      }
      debugTMP = &dynamic_cast<coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>* >(_macroscopicCellServices[0])->getMacroscopicCells();

      preprocessingForMD2Macro(globalCellIndicesFromMacroscopicSolver, size);

      // reset macroscopic data (only those should be used by macroscopic solver anyway) in cells
      for (coupling::datastructures::MacroscopicCell<dim> * & macroscopicCell : _macroscopicCells){
        macroscopicCell->setMacroscopicMass(0.0);
        macroscopicCell->setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
      }

      timeval start{};
      timeval end{};
      double runtime=0;
      if(_rank == 0){ gettimeofday(&start,nullptr); }

      std::vector<coupling::sendrecv::DataExchangeFromMD2Macro<dim> * > allDEs(_totalNumberMDSimulations);
      std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim> *> > allMacroscopicCellsFromMamico(_totalNumberMDSimulations);
      for(unsigned int i=0;i<_totalNumberMDSimulations;++i) {
        allDEs[i] = new coupling::sendrecv::DataExchangeFromMD2Macro<dim> (
          _macroscopicSolverInterface, _macroscopicCellServices[i]->getIndexConversion(),
          _macroscopicCellServices[i]->getID()
        );
      }

      _fromMD2Macro.reduceFromMD2Macro(
        allDEs,
        macroscopicCellsFromMacroscopicSolver,
        globalCellIndicesFromMacroscopicSolver,
        _sumMacroscopicCells
      );

      if(_rank == 0){
        gettimeofday(&end,nullptr);
        runtime += (double)((end.tv_sec-start.tv_sec)*1000000 + (end.tv_usec-start.tv_usec));
      }

      // receive data from each MD simulation and accumulate information in cells
      for (unsigned int l = 0; l < _totalNumberMDSimulations; l++){
        //_macroscopicCellServices[l]->sendFromMD2MacroPostProcess();
      }

      // average data
      for (unsigned int i = 0; i < size; i++) {
        _macroscopicCells[i]->setMacroscopicMass(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMass()/_totalNumberMDSimulations);
        _macroscopicCells[i]->setMacroscopicMomentum(1.0/_totalNumberMDSimulations * macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
      }

      // store data in macroscopicCellsFromMacroscopicSolver
      for (unsigned int i = 0; i < size; i++){
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMass(_macroscopicCells[i]->getMacroscopicMass());
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMomentum(_macroscopicCells[i]->getMacroscopicMomentum());
      }

      debugTMP = &dynamic_cast<coupling::services::MacroscopicCellServiceImpl<LinkedCell,dim>* >(_macroscopicCellServices[0])->getMacroscopicCells();
      if(debugTMP != nullptr) {}
      return runtime;
    }

    /** collects data from MD simulations, averages over them (only macroscopic mass/momentum is considered) and writes the result back into macroscopicCellsFromMacroscopicSolver. */
    double sendFromMD2Macro(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ){
      double res = 0;
      auto size = macroscopicCellsFromMacroscopicSolver.size(); // we assume globalCellIndicesFromMacroscopicSolver to be of identical size

      preprocessingForMD2Macro(globalCellIndicesFromMacroscopicSolver, size);

      // reset macroscopic data (only those should be used by macroscopic solver anyway) in cells
      for (unsigned int i = 0; i < size; i++){
        _macroscopicCells[i]->setMacroscopicMass(0.0);
        _macroscopicCells[i]->setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
        _sumMacroscopicCells[i]->setMacroscopicMass(0.0);
        _sumMacroscopicCells[i]->setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
      }

      // receive data from each MD simulation and accumulate information in cells
      for (unsigned int l = 0; l < _totalNumberMDSimulations; l++){
        //std::cout << "Rank " << _macroscopicCellServices[l]->getIndexConversion().getThisRank() << ": Send from MD to Macro for Simulation no. " << l << std::endl;
        res += _macroscopicCellServices[l]->sendFromMD2Macro(macroscopicCellsFromMacroscopicSolver,globalCellIndicesFromMacroscopicSolver);
        for (unsigned int i = 0; i < size; i++){
          _macroscopicCells[i]->addMacroscopicMass(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMass());
          _macroscopicCells[i]->addMacroscopicMomentum(macroscopicCellsFromMacroscopicSolver[i]->getMacroscopicMomentum());
        }
      }

      // average data
      for (unsigned int i = 0; i < size; i++) {
        _macroscopicCells[i]->setMacroscopicMass(_macroscopicCells[i]->getMacroscopicMass()/_totalNumberMDSimulations);
        _macroscopicCells[i]->setMacroscopicMomentum(1.0/_totalNumberMDSimulations * _macroscopicCells[i]->getMacroscopicMomentum());
      }

      // apply post multi instance FilterPipeline on cell data
      #ifdef ENABLE_POST_MULTI_INSTANCE_FILTERING

      #ifdef DEBUG_FILTER_PIPELINE
      std::cout << "FP: Now applying post-multi-instance filter pipeline" << std::endl;
      #endif

      (*_postMultiInstanceFilterPipeline)();
      #endif  


      // store data in macroscopicCellsFromMacroscopicSolver
      for (unsigned int i = 0; i < size; i++){
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMass(_macroscopicCells[i]->getMacroscopicMass());
        macroscopicCellsFromMacroscopicSolver[i]->setMacroscopicMomentum(_macroscopicCells[i]->getMacroscopicMomentum());
      }

      return res;
    }

    void preprocessingForMD2Macro(
      const unsigned int * const globalCellIndicesFromMacroscopicSolver,
      const unsigned int size
    ) {

      /*
       * If this is first coupling step, we must allocate space for the macroscopic cells we filter and determine averages with.
       */
      if(_cellIndices.empty()) {
        for(unsigned int i = 0; i < size; i++) { _cellIndices.push_back(globalCellIndicesFromMacroscopicSolver[i]); }
      }

      /*
       * If this is the first coupling step, copy cells to be in format compatible with filtering system.
       * Can be optimized if indices don't change dynamically.
       */
      if(_macroscopicCells.empty()) {
        //Allocate & init _macroscopicCells
        for(unsigned int c = _macroscopicCells.size(); c < size; c++) {
          _macroscopicCells.push_back(new coupling::datastructures::MacroscopicCell<dim>());
        }
      }

      /*
       * If this is the first coupling step, we must init the post multi instance filter pipeline operating on averaged cell data.
       * The ENABLE_POST_MULTI_INSTANCE_FILTERING flag is used for debugging purposes and shall be removed later.
       * If you wish to not use post multi-instance filtering in deployment, you can simply leave the corresponding XML-Subtag empty.
       */
#ifdef ENABLE_POST_MULTI_INSTANCE_FILTERING
      if(_postMultiInstanceFilterPipeline == nullptr) {
          //Init filter pipeline
          _postMultiInstanceFilterPipeline = new coupling::FilterPipeline<dim>(_macroscopicCells, _cellIndices, &(_macroscopicCellServices[0]->getIndexConversion()), _macroscopicSolverInterface, _multiMDService, coupling::Scope::postMultiInstance, _filterPipelineConfiguration.c_str());
      }
#endif

      /*
       *  On the first coupling step,  initialize the reduced macroscopicCell vector according to `size`
       */
      if(!_sumMacroscopicCells.empty()) {
        for(unsigned int i = 0 ; i < _sumMacroscopicCells.size(); ++i) {
          _sumMacroscopicCells[i]->setMacroscopicMass(0);
          _sumMacroscopicCells[i]->setMacroscopicMomentum(tarch::la::Vector<dim,double>(0.0));
        }
      }

      /**
       * local reduction step, has to be executed before EVERY sendFromMD2Macro
       */
      sumUpMacroscopicCellsFromMamico();
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
    
    const unsigned int _localNumberMDSimulations; /** number of MD simulations run on the current rank. This can differ for different blocks, i.e. different topologyOffset values. */
    const unsigned int _totalNumberMDSimulations; /** total number of MD simulations */
    coupling::services::MacroscopicCellService<dim> **_macroscopicCellServices; /** pointers of MacroscopicCellService type, one for each MD simulation */
    const unsigned int _topologyOffset; /** topology offset*/
    const int _tws;
    const unsigned int _intNumberProcesses;

    //const coupling::IndexConversion<dim> _indexConversion; /* Used for index conversions during filtering TODO after merge with dynamic-md*/
    std::vector<coupling::datastructures::MacroscopicCell<dim>* > _macroscopicCells; /** used to store in MD data in sendFromMDtoMacro */
    std::vector<coupling::datastructures::MacroscopicCell<dim>* > _sumMacroscopicCells; /** used to reduce all local macroscopic cells before sending from md 2 macro */
    std::vector<unsigned int> _cellIndices; /** used to store in indexing of the above */ 


    /*
     * These are currently only stored in order to properly init the
     * post multi-instance FilterPipeline in sendFromMDToMacro's first invocation.
    */
    coupling::interface::MacroscopicSolverInterface<dim> *_macroscopicSolverInterface;
    const std::string _filterPipelineConfiguration;
    const tarch::utils::MultiMDService<dim>& _multiMDService;
	  const unsigned int _rank;


    /*
     * Analogon to MacroscopicCellService's FilterPipeline. 
     * Is applied during this->sendFromMD2Macro.
    */
    coupling::FilterPipeline<dim> *_postMultiInstanceFilterPipeline;

    coupling::sendrecv::FromMD2Macro<coupling::datastructures::MacroscopicCell<dim>,dim  > _fromMD2Macro;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SERVICES_MULTIMDCELLSERVICE_H_
