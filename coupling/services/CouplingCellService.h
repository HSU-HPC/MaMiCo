// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_COUPLINGCELLSERVICE_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_COUPLINGCELLSERVICE_H_

#pragma once

#include "coupling/CouplingCellPlotter.h"
#include "coupling/KineticEnergyController.h"
#include "coupling/MomentumController.h"
#include "coupling/cell-mappings/ComputeMeanPotentialEnergyMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/PerturbateVelocityMapping.h"
#include "coupling/configurations/BoundaryForceConfiguration.h"
#include "coupling/configurations/CouplingCellConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/configurations/MomentumInsertionConfiguration.h"
#include "coupling/configurations/ParallelTopologyConfiguration.h"
#include "coupling/configurations/ParticleInsertionConfiguration.h"
#include "coupling/configurations/ThermostatConfiguration.h"
#include "coupling/configurations/TransferStrategyConfiguration.h"
#include "coupling/datastructures/CellContainer.h"
#include "coupling/indexing/IndexingService.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/sendrecv/DataExchangeFromMD2Macro.h"
#include "coupling/sendrecv/DataExchangeFromMacro2MD.h"
#include "coupling/sendrecv/FromMD2Macro.h"
#include "coupling/sendrecv/FromMacro2MD.h"
#include "tarch/utils/MultiMDService.h"

#include "coupling/filtering/FilterPipeline.h"

namespace coupling {
namespace services {

template <unsigned int dim> class CouplingCellService;

template <class LinkedCell, unsigned int dim> class CouplingCellServiceImpl;
} // namespace services
} // namespace coupling

/** generic interface class for functionality of data exchange in hybrid
 * Micro-Macro simulations.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::services::CouplingCellService {
public:
  CouplingCellService(unsigned int ID) : _id(ID) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "CouplingCellService::CouplingCellService(): Set ID=" << _id << std::endl;
#endif
  }
  virtual ~CouplingCellService() {}

  virtual double applyFilterPipeline() = 0;
  virtual void sendFromMacro2MDPreProcess() = 0;
  virtual void sendFromMacro2MDPostProcess() = 0;
  virtual void sendFromMacro2MD(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver,
                                const I00* const indices) = 0;
  virtual void sendFromMD2MacroPreProcess() = 0;
  virtual void sendFromMD2MacroPostProcess() = 0;
  virtual double sendFromMD2Macro(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver,
                                  const I00* const indices) = 0;
  virtual void processInnerCouplingCellAfterMDTimestep() = 0;
  virtual void computeAndStoreTemperature(double temperature) = 0;
  virtual void applyTemperatureToMolecules(unsigned int t) = 0;
  virtual void distributeMass(unsigned int t) = 0;
  virtual void distributeMomentum(unsigned int t) = 0;
  virtual void applyBoundaryForce(unsigned int t) = 0;
  virtual void perturbateVelocity() = 0;
  virtual void plotEveryMicroscopicTimestep(unsigned int t) = 0;
  virtual void plotEveryMacroscopicTimestep(unsigned int t) = 0;

  virtual void initFiltering() {
    throw std::runtime_error("CouplingCellService: Error: Called "
                             "initFiltering for non-Impl object.");
  } /*Note: This is not pure virtual, because some implementations of this
       interface don't have a FilterPipeline. */
  virtual const coupling::filtering::FilterPipeline<dim>* getFilterPipeline() const {
    throw std::runtime_error("CouplingCellService: Error: Called getFilterPipeline() in instance "
                             "without FilterPipeline.");
  } /*Note: This is not pure virtual, because some implementations of this
       interface don't have a FilterPipeline. */

  unsigned int getID() const { return _id; }

protected:
  const unsigned int _id; /** (unique) identifier of this coupling cell service */
};

/** This class put together all ingredients for coupling MD and some macroscopic
 * solver. It thus triggers send/recv-operations between the coupling tool and
 * MD as well as between the coupling tool and the macroscopic solver.
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::services::CouplingCellServiceImpl : public coupling::services::CouplingCellService<dim> {
public:
  /** constructor. Arguments:
   *  mdSolverInterface                 - pointer to the MD solver interface
   *  macroscopicSolverInterface        - pointer to the macroscopic solver
   * interface numberProcesses                   - number of processes assuming
   * a block-like domain decomposition rank                              - rank
   * of the current process particleInsertionConfiguration    - configuration
   * object for (USHER-based) particle insertion and particle removal
   *  momentumInsertionConfiguration    - configuration object which determines
   * the momentum transfer on MD side transferStrategyConfiguration     -
   * configuration object which determines the respective transfer strategy
   *  parallelTopologyConfiguration     - configuratio object which defines the
   * parallel topology of the simulation (domain decomposition of the MD
   * simulation) numberMDTimestepsPerCouplingCycle - number of MD time steps per
   * coupling cycle couplingCellConfiguration      - configuration object
   * which determines the properties of the coupling cells topologyOffset -
   * offset in linearized topology of ranks topologyGlobalNumberProcesses     -
   * global number of processes available in overall topology
   *
   *  Note: the interface pointers are used by the CouplingCellServiceImpl;
   * they are not deleted at the end of the simulation since other (external)
   * routines may use those as well.
   */
  CouplingCellServiceImpl(
      unsigned int ID, coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface,       // interface to MD simulation
      coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface,                  // interface to macroscopic solver
      tarch::la::Vector<dim, unsigned int> numberProcesses,                                              // number of processes in all directions
      unsigned int rank,                                                                                 // current rank
      const coupling::configurations::ParticleInsertionConfiguration& particleInsertionConfiguration,    // configuration for particle
                                                                                                         // insertion
      const coupling::configurations::MomentumInsertionConfiguration& momentumInsertionConfiguration,    // configuration for momentum
                                                                                                         // insertion
      const coupling::configurations::BoundaryForceConfiguration<dim>& boundaryForceConfiguration,       // configuration for boundary forces
      const coupling::configurations::TransferStrategyConfiguration<dim>& transferStrategyConfiguration, // configuration for transfer strategy
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,      // configuration for parallel topology
      const coupling::configurations::ThermostatConfiguration& thermostatConfiguration,
      unsigned int numberMDTimestepsPerCouplingCycle,                                            // number MD timesteps per coupling
                                                                                                 // cycle (required to initialise
                                                                                                 // transfer strategy)
      const coupling::configurations::CouplingCellConfiguration<dim>& couplingCellConfiguration, // configuration for coupling cells
                                                                                                 // and respective plotting
      const char* filterPipelineConfiguration, const tarch::utils::MultiMDService<dim>& multiMDService, unsigned int topologyOffset, int tws = 0);

  CouplingCellServiceImpl(
      unsigned int ID, coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface,       // interface to MD simulation
      coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface,                  // interface to macroscopic solver
      tarch::la::Vector<dim, unsigned int> numberProcesses,                                              // number of processes in all directions
      unsigned int rank,                                                                                 // current rank
      const coupling::configurations::ParticleInsertionConfiguration& particleInsertionConfiguration,    // configuration for particle
                                                                                                         // insertion
      const coupling::configurations::MomentumInsertionConfiguration& momentumInsertionConfiguration,    // configuration for momentum
                                                                                                         // insertion
      const coupling::configurations::BoundaryForceConfiguration<dim>& boundaryForceConfiguration,       // configuration for boundary forces
      const coupling::configurations::TransferStrategyConfiguration<dim>& transferStrategyConfiguration, // configuration for transfer strategy
      const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration,      // configuration for parallel topology
      const coupling::configurations::ThermostatConfiguration& thermostatConfiguration,
      unsigned int numberMDTimestepsPerCouplingCycle,                                            // number MD timesteps per coupling
                                                                                                 // cycle (required to initialise
                                                                                                 // transfer strategy)
      const coupling::configurations::CouplingCellConfiguration<dim>& couplingCellConfiguration, // configuration for coupling cells
                                                                                                 // and respective plotting
      const char* filterPipelineConfiguration, const tarch::utils::MultiMDService<dim>& multiMDService)
      : CouplingCellServiceImpl<LinkedCell, dim>(ID, mdSolverInterface, macroscopicSolverInterface, numberProcesses, rank, particleInsertionConfiguration,
                                                 momentumInsertionConfiguration, boundaryForceConfiguration, transferStrategyConfiguration,
                                                 parallelTopologyConfiguration, thermostatConfiguration, numberMDTimestepsPerCouplingCycle,
                                                 couplingCellConfiguration, filterPipelineConfiguration, multiMDService, 0) {}

  /** destructor. Frees dynamically allocated memory for particle insertion,
   * momentum insertion and the transfer strategy.
   */
  virtual ~CouplingCellServiceImpl();

  /** sends information from macroscopic solver to MD. The cell information from
   * the macroscopic solver is handed over as array including corresponding
   * global cell indices and the number of cells. The coupling tool internally
   * triggers send/recv-operations (this also comprises the distributed memory
   *  parallelisation via MPI) and writes the respective information to the
   * coupling cells of the tool
   */
  void sendFromMacro2MD(const coupling::datastructures::FlexibleCellContainer<dim>& macroscopicSolverCellContainer) override;

  /** sends information from MD to the macroscopic solver. After the
   * send/recv-operations (this also comprises the distributed memory
   * parallelisation scenario), the information from the coupling tool is
   * written to the buffer couplingCellsFromMacroscopicSolver together with
   * the respective global cell indices
   *  (-> indices).
   * 
   * @returns The runtime of filtering related code in usec.
   */
  double sendFromMD2Macro(const coupling::datastructures::FlexibleCellContainer& macroscopicSolverCellContainer) override;

  /** applies the filter pipeline and returns the runtime of this operation */
  double applyFilterPipeline() override;

  void sendFromMacro2MDPreProcess() override;

  void sendFromMacro2MDPostProcess() override;

  void sendFromMD2MacroPreProcess() override;

  void sendFromMD2MacroPostProcess() override {}

  /** carries out coupling-dependent operations (such as sampling) on the
   * non-ghost coupling cells after each MD time step. This method needs thus
   * to be called from the MD simulation.
   */
  void processInnerCouplingCellAfterMDTimestep() override;

  /** sets the temperature value in all coupling cells. If the value of
   * temperature is -1.0, we set the local temperature of each coupling cell
   * (and just store this value in the coupling cell). Otherwise, we apply
   * the given temperature in all cells. In the latter case, this also resembles
   * a first thermostat-like operation.
   */
  void computeAndStoreTemperature(double temperature) override;

  /** applies a thermostat in all non-ghost coupling cells. */
  void applyTemperatureToMolecules(unsigned int t) override;

  /** applies a boundary force to molecules which are close to an open boundary.
   */
  void applyBoundaryForce(unsigned int t) override;

  /** distributes mass in the system. */
  void distributeMass(unsigned int t) override;

  /** distributes momentum in MD. Should typically be called after force
   * accumulation since momentum distribution may depend on current forces. */
  void distributeMomentum(unsigned int t) override;

  /** applies a new velocity to each particle according to its cell's mean
   * velocity. */
  void perturbateVelocity() override;

  /** plots coupling cell and molecule information at some time step t. The
   * correct triggering of plotting needs to be established from the main
   * coupling loop which is outside the coupling tool (not included in this
   * function).
   */
  void plotEveryMicroscopicTimestep(unsigned int t) override;
  void plotEveryMacroscopicTimestep(unsigned int t) override;

  /**
   * Initialises the _filterPipeline member. Called from _multiMDCellService's
   * constructFilterPipelines(). Make sure to delete _filterPipeline in
   * ~CouplingCellServiceImpl()
   */
  void initFiltering() override {
    _filterPipeline = new coupling::filtering::FilterPipeline<dim>(_couplingCells.getCouplingCells(), coupling::filtering::Scope::perInstance, _multiMDService,
                                                                   _filterPipelineConfiguration);
  }

  const coupling::filtering::FilterPipeline<dim>* getFilterPipeline() const override { return _filterPipeline; }

  /**
   * Creates a new filter from scratch and appends it to a sequence that is part
   * of this service's filter pipelining system. For that, the desired
   * sequence's identifier and two functions are needed:
   *  - applyScalar What to do with scalar properties of the sequence's
   * Coupling Cells.
   *  - applyVector: What to do with properties stored as vectors of the
   * sequence's of Coupling Cells.
   */
  /*
   * TODO: MOVE COMMENT
  void addFilterToSequence(	const char *sequenceIdentifier,
                                                          const
  std::function<std::vector<double> (std::vector<double> cells_s,
  std::vector<std::array<unsigned int, dim>> indices)>* applyScalar, const
  std::function<std::vector<std::array<double, dim>>
  (std::vector<std::array<double, dim>> cells_v, std::vector<std::array<unsigned
  int, dim>> indices)>* applyVector, int filterIndex = -1
  );*/

  /** returns the coupling cells. This functions is meant to be used in test
   * scenarios and for debugging only! DO NOT USE IT FOR OTHER PURPOSES! */
  coupling::datastructures::LinkedCellContainer<LinkedCell, dim>& getCouplingCells() { return _couplingCells; }

private:
// ------------------- INCLUDE WRAPPER DEFINITIONS
// -------------------------------------
#include "CouplingCellTraversalWrappers.cpph"

  /** initialises the index structures for USHER scheme */
  void initIndexVectors4Usher(tarch::la::Vector<dim, unsigned int> numberLinkedCellsPerCouplingCell);

  tarch::la::Vector<dim, double> getPositionOfFirstLocalGhostCell() const;

  std::function<void(Wrapper&)> initCorrectApplicationOfThermostat(const coupling::configurations::ThermostatConfiguration& thermostatConfiguration) {
    if (thermostatConfiguration.getThermostatRegionType() == coupling::configurations::ThermostatConfiguration::ThermostatRegion::all)
      return [this](Wrapper& wrapper) { _couplingCells.applyToLocalNonGhostCouplingCellsWithLinkedCells(wrapper); };
    else if (thermostatConfiguration.getThermostatRegionType() == coupling::configurations::ThermostatConfiguration::ThermostatRegion::outerLayers)
      return [this, &thermostatConfiguration](Wrapper& wrapper) {
        _couplingCells.applyXLayersOfGlobalNonGhostCellsWithLinkedCells(wrapper, thermostatConfiguration.getCells2Use());
      };
    else if (thermostatConfiguration.getThermostatRegionType() == coupling::configurations::ThermostatConfiguration::ThermostatRegion::onlyOutestLayer)
      return [this](Wrapper& wrapper) { _couplingCells.applyToFirstLayerOfGlobalNonGhostCellsWithLinkedCells(wrapper); };
    else
      return [](Wrapper& wrapper) {};
  }

  /** number of MD time steps in each coupling cycle */
  const unsigned int _numberMDTimestepsPerCouplingCycle;

  /** interface for MD solver */
  coupling::interface::MDSolverInterface<LinkedCell, dim>* _mdSolverInterface;

  /** interface for macroscopic solver */
  coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface;

  /** for quantity transfer between solvers */
  coupling::sendrecv::FromMacro2MD<coupling::datastructures::CouplingCell<dim>, dim> _fromMacro2MD;
  coupling::sendrecv::DataExchangeFromMacro2MD<dim> _deFromMacro2MD;
  coupling::sendrecv::FromMD2Macro<coupling::datastructures::CouplingCell<dim>, dim> _fromMD2Macro;
  coupling::sendrecv::DataExchangeFromMD2Macro<dim> _deFromMD2Macro;

  /** storage for coupling cells in coupling tool */
  coupling::datastructures::LinkedCellContainer<LinkedCell, dim> _couplingCells;

  /** filter pipeline, used to apply filters in sendFromMD2Macro */
  coupling::filtering::FilterPipeline<dim>* _filterPipeline;

  /**parameters needed in initFiltering() */
  const char* _filterPipelineConfiguration;
  const tarch::utils::MultiMDService<dim> _multiMDService;

  /** needed for insertion of momentum */
  coupling::MomentumInsertion<LinkedCell, dim>* _momentumInsertion;
  coupling::configurations::MomentumInsertionConfiguration::MomentumInsertionType _momentumInsertionType;
  /** needed for insertion of particles, e.g. USHER */
  coupling::ParticleInsertion<LinkedCell, dim>* _particleInsertion;
  const tarch::la::Vector<dim, unsigned int> _numberLinkedCellsPerCouplingCell;
  const coupling::configurations::ParticleInsertionConfiguration::ParticleInsertionType _particleInsertionType;
  /** coupling strategy */
  coupling::transferstrategies::TransferStrategy<LinkedCell, dim>* _transferStrategy;

  /** controls the kinetic energy of the system, i.e. maintains temperature in
   * case of changing mass/momentum. */
  coupling::KineticEnergyController<LinkedCell, dim> _kineticEnergyController;
  /** controls and apply boundary forces to molecules close to open
   * boundaries.*/
  coupling::BoundaryForceController<LinkedCell, dim>* _boundaryForceController;
  /** controls/ maintains momentum, e.g. after particle insertion */
  coupling::MomentumController<LinkedCell, dim> _momentumController;

  std::function<void(Wrapper&)> _applyAccordingToConfiguration;

  /** information for plotting */
  const std::string _microscopicFilename;
  const unsigned int _writeEveryMicroscopicTimestep;
  const std::string _macroscopicFilename;
  const unsigned int _writeEveryMacroscopicTimestep;

  /** index vectors for block-usher scheme
   * -----------------------------------------------------*/
  // start and end coordinate for block loop over coupling cells (with 3
  // entries always!)
  tarch::la::Vector<3, unsigned int> _usherCellStart[1 << dim];
  tarch::la::Vector<3, unsigned int> _usherCellEnd[1 << dim];
  tarch::la::Vector<dim, unsigned int> _usherRange[1 << dim];
  // offset in red-black loops nested within the block loops (always 0 or 1
  // entries)
  tarch::la::Vector<3, unsigned int> _usherCellOffset[1 << dim];
};
#include "CouplingCellService.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_SERVICES_COUPLINGCELLSERVICE_H_
