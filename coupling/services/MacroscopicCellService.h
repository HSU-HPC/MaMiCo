// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_MACROSCOPICCELLSERVICE_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_MACROSCOPICCELLSERVICE_H_

#pragma once

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/sendrecv/FromMacro2MD.h"
#include "coupling/sendrecv/DataExchangeFromMacro2MD.h"
#include "coupling/sendrecv/FromMD2Macro.h"
#include "coupling/sendrecv/DataExchangeFromMD2Macro.h"
#include "coupling/datastructures/MacroscopicCells.h"
#include "coupling/MacroscopicCellPlotter.h"
#include "coupling/configurations/ParticleInsertionConfiguration.h"
#include "coupling/configurations/MomentumInsertionConfiguration.h"
#include "coupling/configurations/BoundaryForceConfiguration.h"
#include "coupling/configurations/TransferStrategyConfiguration.h"
#include "coupling/configurations/ParallelTopologyConfiguration.h"
#include "coupling/configurations/MacroscopicCellConfiguration.h"
#include "coupling/KineticEnergyController.h"
#include "coupling/MomentumController.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/ComputeMeanPotentialEnergyMapping.h"
#include "tarch/utils/MultiMDService.h"

#include "coupling/filtering/FilterPipeline.h"

namespace coupling {
  namespace services {

    template<unsigned int dim>
    class MacroscopicCellService;

    template<class LinkedCell,unsigned int dim>
    class MacroscopicCellServiceImpl;
  }
}



/** generic interface class for functionality of data exchange in hybrid Micro-Macro simulations.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::services::MacroscopicCellService {
  public:
    MacroscopicCellService(unsigned int ID): _id(ID){
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "MacroscopicCellService::MacroscopicCellService(): Set ID=" << _id << std::endl;
      #endif
    }
    virtual ~MacroscopicCellService(){}

    virtual void sendFromMacro2MD(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ) = 0;
    virtual double sendFromMD2Macro(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    ) = 0;
    virtual void sendFromMacro2MDPreProcess() = 0;
    virtual void sendFromMacro2MDPostProcess() = 0;
    virtual void processInnerMacroscopicCellAfterMDTimestep() = 0;
    virtual void computeAndStoreTemperature(double temperature) = 0;
    virtual void applyTemperatureToMolecules(unsigned int t) = 0;
    virtual void distributeMass(unsigned int t) = 0;
    virtual void distributeMomentum(unsigned int t) = 0;
    virtual void applyBoundaryForce(unsigned int t) = 0;
    virtual void plotEveryMicroscopicTimestep(unsigned int t) = 0;
    virtual void plotEveryMacroscopicTimestep(unsigned int t) = 0;
    virtual const coupling::IndexConversion<dim>& getIndexConversion() const = 0;
	virtual const coupling::FilterPipeline<dim>& getFilterPipeline() const { throw std::runtime_error("MacroscopicCellService: Error: Called getFilterPipeline() in instance without FilterPipeline."); }  /*Note: This is not pure virtual, because some implementations of this interface don't have a FilterPipeline. */
	unsigned int getID() const { return _id;}


  protected:
    const unsigned int _id; /** (unique) identifier of this macroscopic cell service */
};

/** This class put together all ingredients for coupling MD and some macroscopic solver.
 *  It thus triggers send/recv-operations between the coupling tool and MD as well as between the coupling tool
 *  and the macroscopic solver.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::services::MacroscopicCellServiceImpl:
public coupling::services::MacroscopicCellService<dim> {
  public:
    /** constructor. Arguments:
     *  mdSolverInterface                 - pointer to the MD solver interface
     *  macroscopicSolverInterface        - pointer to the macroscopic solver interface
     *  numberProcesses                   - number of processes assuming a block-like domain decomposition
     *  rank                              - rank of the current process
     *  particleInsertionConfiguration    - configuration object for (USHER-based) particle insertion and particle removal
     *  momentumInsertionConfiguration    - configuration object which determines the momentum transfer on MD side
     *  transferStrategyConfiguration     - configuration object which determines the respective transfer strategy
     *  parallelTopologyConfiguration     - configuratio object which defines the parallel topology of the simulation (domain decomposition of the MD simulation)
     *  numberMDTimestepsPerCouplingCycle - number of MD time steps per coupling cycle
     *  macroscopicCellConfiguration      - configuration object which determines the properties of the macroscopic cells
     *  topologyOffset                    - offset in linearized topology of ranks
     *  topologyGlobalNumberProcesses     - global number of processes available in overall topology
     *
     *  Note: the interface pointers are used by the MacroscopicCellServiceImpl; they are not deleted at the end of the
     *  simulation since other (external) routines may use those as well.
     */
    MacroscopicCellServiceImpl(
      unsigned int ID,
      coupling::interface::MDSolverInterface<LinkedCell,dim> *mdSolverInterface,                       // interface to MD simulation
      coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface,                // interface to macroscopic solver
      tarch::la::Vector<dim,unsigned int> numberProcesses,                                                                // number of processes in all directions
      unsigned int rank,                                                                                                  // current rank
      const coupling::configurations::ParticleInsertionConfiguration &particleInsertionConfiguration,  // configuration for particle insertion
      const coupling::configurations::MomentumInsertionConfiguration &momentumInsertionConfiguration,  // configuration for momentum insertion
      const coupling::configurations::BoundaryForceConfiguration<dim> &boundaryForceConfiguration,     // configuration for boundary forces
      const coupling::configurations::TransferStrategyConfiguration<dim>& transferStrategyConfiguration,    // configuration for transfer strategy
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,    // configuration for parallel topology
      unsigned int numberMDTimestepsPerCouplingCycle,                                                  // number MD timesteps per coupling cycle (required to initialise transfer strategy)
      const coupling::configurations::MacroscopicCellConfiguration<dim> &macroscopicCellConfiguration, // configuration for macroscopic cells and respective plotting
	  const char* filterPipelineConfiguration,
      const tarch::utils::MultiMDService<dim>& multiMDService,
      unsigned int topologyOffset, int tws = 0
    );

    MacroscopicCellServiceImpl(
      unsigned int ID,
      coupling::interface::MDSolverInterface<LinkedCell,dim> *mdSolverInterface,                       // interface to MD simulation
      coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface,                // interface to macroscopic solver
      tarch::la::Vector<dim,unsigned int> numberProcesses,                                                                // number of processes in all directions
      unsigned int rank,                                                                                                  // current rank
      const coupling::configurations::ParticleInsertionConfiguration &particleInsertionConfiguration,  // configuration for particle insertion
      const coupling::configurations::MomentumInsertionConfiguration &momentumInsertionConfiguration,  // configuration for momentum insertion
      const coupling::configurations::BoundaryForceConfiguration<dim> &boundaryForceConfiguration,     // configuration for boundary forces
      const coupling::configurations::TransferStrategyConfiguration<dim>& transferStrategyConfiguration,    // configuration for transfer strategy
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,    // configuration for parallel topology
      unsigned int numberMDTimestepsPerCouplingCycle,                                                  // number MD timesteps per coupling cycle (required to initialise transfer strategy)
      const coupling::configurations::MacroscopicCellConfiguration<dim> &macroscopicCellConfiguration,  // configuration for macroscopic cells and respective plotting
	  const char* filterPipelineConfiguration,
      const tarch::utils::MultiMDService<dim>& multiMDService
    ): MacroscopicCellServiceImpl<LinkedCell,dim>(ID,mdSolverInterface,macroscopicSolverInterface,numberProcesses,rank,
       particleInsertionConfiguration,momentumInsertionConfiguration,boundaryForceConfiguration,transferStrategyConfiguration,
       parallelTopologyConfiguration,numberMDTimestepsPerCouplingCycle,macroscopicCellConfiguration,multiMDService,0){}

    /** destructor. Frees dynamically allocated memory for particle insertion, momentum insertion and the transfer
     *  strategy.
     */
    virtual ~MacroscopicCellServiceImpl();

    /** sends information from macroscopic solver to MD. The cell information from the macroscopic solver is
     *  handed over as array including corresponding global cell indices and the number of cells.
     *  The coupling tool internally triggers send/recv-operations (this also comprises the distributed memory
     *  parallelisation via MPI) and writes the respective information to the macroscopic cells of the tool
     */
    void sendFromMacro2MD(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    );


    /** sends information from MD to the macroscopic solver. After the send/recv-operations (this also comprises
     *  the distributed memory parallelisation scenario), the information from the coupling tool is written to the
     *  buffer macroscopicCellsFromMacroscopicSolver together with the respective global cell indices
     *  (-> globalCellIndicesFromMacroscopicSolver).
     */
    double sendFromMD2Macro(
      const std::vector<coupling::datastructures::MacroscopicCell<dim>* > &macroscopicCellsFromMacroscopicSolver,
      const unsigned int * const globalCellIndicesFromMacroscopicSolver
    );

    void sendFromMacro2MDPreProcess();

    void sendFromMacro2MDPostProcess();


    /** carries out coupling-dependent operations (such as sampling) on the non-ghost macroscopic cells after each
     *  MD time step. This method needs thus to be called from the MD simulation.
     */
    void processInnerMacroscopicCellAfterMDTimestep();

    /** sets the temperature value in all macroscopic cells. If the value of temperature is -1.0,
     *  we set the local temperature of each macroscopic cell (and just store this value in the macroscopic cell).
     *  Otherwise, we apply the given temperature in all cells. In the latter case, this also resembles a first
     *  thermostat-like operation.
     */
    void computeAndStoreTemperature(double temperature);


    /** applies a thermostat in all non-ghost macroscopic cells. */
    void applyTemperatureToMolecules(unsigned int t);

    /** applies a boundary force to molecules which are close to an open boundary. */
    void applyBoundaryForce(unsigned int t);

    /** distributes mass in the system. */
    void distributeMass(unsigned int t);

    /** distributes momentum in MD. Should typically be called after force accumulation since momentum distribution may depend on current forces. */
    void distributeMomentum(unsigned int t);

    /** plots macroscopic cell and molecule information at some time step t. The correct triggering of plotting needs
     *  to be established from the main coupling loop which is outside the coupling tool (not included in this function).
     */
    void plotEveryMicroscopicTimestep(unsigned int t);
    void plotEveryMacroscopicTimestep(unsigned int t);

    /** returns a reference to the index conversion object. Some external classes may require information on macroscopic
     *  cell size etc. which can then be easily accessed. This is thus not required by the internal mechanisms of the
     *  coupling tool.
     */
    const coupling::IndexConversion<dim>& getIndexConversion() const { return *_indexConversion; }
	
	const coupling::FilterPipeline<dim>& getFilterPipeline() const { return _filterPipeline; }

	/**
	 * Creates a new filter from scratch and appends it to a sequence that is part of this service's filter pipelining system.
	 * For that, the desired sequence's identifier and two functions are needed:
	 *  - applyScalar What to do with scalar properties of the sequence's Macroscopic Cells.
	 *  - applyVector: What to do with properties stored as vectors of the sequence's of Macroscopic Cells.
	 */
	/*
	 * TODO: MOVE COMMENT
	void addFilterToSequence(	const char *sequenceIdentifier,
		   						const std::function<std::vector<double> (std::vector<double> cells_s, std::vector<std::array<unsigned int, dim>> indices)>* applyScalar,
		   						const std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>> cells_v, std::vector<std::array<unsigned int, dim>> indices)>* applyVector,
								int filterIndex = -1
	);*/
	
    /** returns the macroscopic cells. This functions is meant to be used in test scenarios and for debugging only! DO NOT USE IT FOR OTHER PURPOSES! */
    coupling::datastructures::MacroscopicCells<LinkedCell,dim>& getMacroscopicCells() { return _macroscopicCells;}
	
  private:
    /** initialises the IndexConversion object at start up. This is the very first thing to be done in the
     *  constructor since nearly all subsequent operations depend on indexing of cells.
     */
    coupling::IndexConversion<dim>* initIndexConversion(
    tarch::la::Vector<dim,double> macroscopicCellSize, tarch::la::Vector<dim,unsigned int>numberProcesses, unsigned int rank,
    tarch::la::Vector<dim,double> globalMDDomainSize, tarch::la::Vector<dim,double> globalMDDomainOffset,
    coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
    unsigned int topologyOffset) const;

    /** initialises the index structures for USHER scheme */
    void initIndexVectors4Usher( tarch::la::Vector<dim,unsigned int> numberLinkedCellsPerMacroscopicCell );
	
    tarch::la::Vector<dim,double> getPositionOfFirstLocalGhostCell() const;
	


    /** needed to determine cell range, ranks etc. */
    const coupling::IndexConversion<dim> *_indexConversion;
    /** number of MD time steps in each coupling cycle */
    const unsigned int _numberMDTimestepsPerCouplingCycle;

    /** interface for MD solver */
    coupling::interface::MDSolverInterface<LinkedCell,dim> *_mdSolverInterface;

    /** interface for macroscopic solver */
    coupling::interface::MacroscopicSolverInterface<dim> *_macroscopicSolverInterface;

    /** for quantity transfer between solvers */
    coupling::sendrecv::FromMacro2MD<coupling::datastructures::MacroscopicCell<dim>,dim > _fromMacro2MD;
    coupling::sendrecv::DataExchangeFromMacro2MD<dim> _deFromMacro2MD;
    coupling::sendrecv::FromMD2Macro<coupling::datastructures::MacroscopicCell<dim>,dim > _fromMD2Macro;
    coupling::sendrecv::DataExchangeFromMD2Macro<dim> _deFromMD2Macro;

    /** storage for macroscopic cells in coupling tool */
    coupling::datastructures::MacroscopicCells<LinkedCell,dim> _macroscopicCells;

    /** filter pipeline, used to apply filters in sendFromMD2Macro */
    coupling::FilterPipeline<dim> _filterPipeline;

    /** needed for insertion of momentum */
    coupling::MomentumInsertion<LinkedCell,dim>* _momentumInsertion;
    coupling::configurations::MomentumInsertionConfiguration::MomentumInsertionType _momentumInsertionType;
    /** needed for insertion of particles, e.g. USHER */
    coupling::ParticleInsertion<LinkedCell,dim>* _particleInsertion;
    const tarch::la::Vector<dim,unsigned int> _numberLinkedCellsPerMacroscopicCell;
    const coupling::configurations::ParticleInsertionConfiguration::ParticleInsertionType _particleInsertionType;
    /** coupling strategy */
    coupling::transferstrategies::TransferStrategy<LinkedCell,dim>* _transferStrategy;

    /** controls the kinetic energy of the system, i.e. maintains temperature in case of changing mass/momentum. */
    coupling::KineticEnergyController<LinkedCell,dim> _kineticEnergyController;
    /** controls and apply boundary forces to molecules close to open boundaries.*/
    coupling::BoundaryForceController<LinkedCell,dim> *_boundaryForceController;
    /** controls/ maintains momentum, e.g. after particle insertion */
    coupling::MomentumController<LinkedCell,dim> _momentumController;

    /** information for plotting */
    const std::string _microscopicFilename;
    const unsigned int _writeEveryMicroscopicTimestep;
    const std::string _macroscopicFilename;
    const unsigned int _writeEveryMacroscopicTimestep;

	//Inner cells managed by std::vectors. Indexing starts at the bottom left inner cell.
	//TODO: REMOVE std::vector<coupling::datastructures::MacroscopicCell<dim> *> _innerMacroscopicCells;
	//TODO: REMOVE std::vector<tarch::la::Vector<dim, unsigned int>> _innerMacroscopicCellIndices;

    /** index vectors for block-usher scheme -----------------------------------------------------*/
    // start and end coordinate for block loop over macroscopic cells (with 3 entries always!)
    tarch::la::Vector<3,unsigned int> _usherCellStart[1<<dim];
    tarch::la::Vector<3,unsigned int> _usherCellEnd[1<<dim];
    tarch::la::Vector<dim,unsigned int> _usherRange[1<<dim];
    // offset in red-black loops nested within the block loops (always 0 or 1 entries)
    tarch::la::Vector<3,unsigned int> _usherCellOffset[1<<dim];

    // ------------------- INCLUDE WRAPPER DEFINITIONS -------------------------------------
    #include "MacroscopicCellTraversalWrappers.cpph"
};
#include "MacroscopicCellService.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_MACROSCOPICCELLSERVICE_H_
