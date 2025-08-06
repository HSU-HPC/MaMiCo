#pragma once
#include "coupling/interface/MDSimulation.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "simplemd/MolecularDynamicsSimulation.h"
#include "simplemd/molecule-mappings/WriteCheckPointMapping.h"

namespace coupling {
namespace interface {
/** defines a MDsimulation based on simplemd::MolecularDynamicsSimulation but for coupled MD simulations.
 * Thus, the implementation of one timestep slightly differs from the one of the base class.
 *	@brief defines MD simulation from default simple MD code.
 *  @author Philipp Neumann
 */
class SimpleMDSimulation : public coupling::interface::MDSimulation, public simplemd::MolecularDynamicsSimulation {
public:
  SimpleMDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration)
      : coupling::interface::MDSimulation(), simplemd::MolecularDynamicsSimulation(configuration), _couplingCellService(NULL), _couplingSwitchedOn(true) {}

  virtual ~SimpleMDSimulation() {}

  virtual void switchOffCoupling() { _couplingSwitchedOn = false; }

  virtual void switchOnCoupling() { _couplingSwitchedOn = true; }

  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) {
    for (unsigned int t = firstTimestep; t < firstTimestep + numberTimesteps; t++) {
      simulateOneCouplingTimestep(t);
    }
  }

  void simulateOneCouplingTimestep(const unsigned int& t) {
    // if coupling is switched off, perform "normal" MD timestep
    if (!_couplingSwitchedOn) {
      simulateOneTimestep(t);
      return;
    }
    if (_parallelTopologyService->isIdle()) {
      return;
    }

    _boundaryTreatment->putBoundaryParticlesToInnerCellsAndFillBoundaryCells(_localBoundary, *_parallelTopologyService);

    // call to synchronise data in cells; needs to be at this point of the
    // coupling algorithm as the particles need to be placed inside the correct
    // sampling volumes (hence: after communication with neighbours and molecule
    // updates); do it BEFORE quantities are manipulated as we can then also do
    // some pre-processing here.
    _couplingCellService->processInnerCouplingCellAfterMDTimestep();
    // ------------ coupling step: distribute mass ---------------------
    _couplingCellService->distributeMass(t);
    // for isothermal simulations: apply thermostat
    _couplingCellService->applyTemperatureToMolecules(t);

    // ---------- from here: go on with usual MD algorithm
    // ------------------------------

    // compute forces. After this step, each molecule has received all force
    // contributions from its neighbors.
    _linkedCellService->iterateCellPairs(*_lennardJonesForce);

    // distribute momentum -> some methods require modification of force terms,
    // therefore we call it AFTER the force computation and before everything else
    _couplingCellService->distributeMomentum(t);
    // apply boundary forces
    _couplingCellService->applyBoundaryForce(t);
    // evaluate statistics
    evaluateStatistics(t);

    _boundaryTreatment->emptyGhostBoundaryCells();

    // plot VTK output
    if ((_configuration.getVTKConfiguration().getWriteEveryTimestep() != 0) && (t % _configuration.getVTKConfiguration().getWriteEveryTimestep() == 0)) {
      _vtkMoleculeWriter->setTimestep(t);
      _moleculeService->getContainer().iterateMolecules(*_vtkMoleculeWriter);
    }

// plot ADIOS2 output
#if BUILD_WITH_ADIOS2
    if ((_configuration.getAdios2Configuration().getWriteEveryTimestep() != 0) && (t % _configuration.getAdios2Configuration().getWriteEveryTimestep() == 0)) {
      _Adios2Writer->setTimestep(t);
      _moleculeService->getContainer().iterateMolecules(*_Adios2Writer);
    }
#endif

    // write checkpoint
    if ((_configuration.getCheckpointConfiguration().getWriteEveryTimestep() != 0) &&
        (t % _configuration.getCheckpointConfiguration().getWriteEveryTimestep() == 0)) {
      simplemd::moleculemappings::WriteCheckPointMapping writeCheckPointMapping(*_parallelTopologyService, _configuration.getCheckpointConfiguration().getFilename(), t);
      _moleculeService->getContainer().iterateMolecules(writeCheckPointMapping);
    }
    // reorganise memory if needed
    if ((_configuration.getSimulationConfiguration().getReorganiseMemoryEveryTimestep() != 0) &&
        (t % _configuration.getSimulationConfiguration().getReorganiseMemoryEveryTimestep() == 0)) {
      _moleculeService->getContainer().sort();
    }
    // plot also coupling cell information
    _couplingCellService->plotEveryMicroscopicTimestep(t);

    // time integration. After this step, the velocities and the positions of the
    // molecules have been updated.
    _moleculeService->getContainer().iterateMolecules(*_timeIntegrator);

    // sort molecules into linked cells
    _moleculeService->getContainer().sort();

    if (_parallelTopologyService->getProcessCoordinates() == tarch::la::Vector<MD_DIM, unsigned int>(0)) {
      // if(t%50==0) std::cout <<"Finish MD timestep " << t << "..." << std::endl;
    }
  }

  virtual void sortMoleculesIntoCells() {
    // nop required, since the linked cells are very tightly linked to mamico
  }

  virtual void setCouplingCellService(coupling::services::CouplingCellService<MDSIMULATIONFACTORY_DIMENSION>* couplingCellService) {
    _couplingCellService = couplingCellService;
    // set the cell service also in singleton of mamico interface provider ->
    // typically not required in coupling, but makes the simulation state more
    // consistent compared to using LAMMPS
    coupling::interface::MamicoInterfaceProvider<simplemd::LinkedCell, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setCouplingCellService(
        couplingCellService);
  }

  virtual void init() { initServices(); }

  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) {
    initServices(multiMDService, localMDSimulation);
  }

  virtual void shutdown() { shutdownServices(); }

  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) {
    simplemd::moleculemappings::WriteCheckPointMapping writeCheckPointMapping(getParallelTopologyService(), filestem, t);
    _moleculeService->getContainer().iterateMolecules(writeCheckPointMapping);
  }

  // function particularly needed to init MD solver interface -> should only be
  // called from factory
  simplemd::BoundaryTreatment& getBoundaryTreatment() { return *_boundaryTreatment; }
  simplemd::services::ParallelTopologyService& getParallelTopologyService() { return *_parallelTopologyService; }
  simplemd::services::MoleculeService& getMoleculeService() {
    #if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if(_moleculeService == NULL){
      std::cout <<"ERROR coupling::interface::MDSimulation::getMoleculeService(): _moleculeService == NULL " << std::endl;
      exit(1);
    }
    #endif
    return *_moleculeService;
  }
  simplemd::services::LinkedCellService& getLinkedCellService() { return *_linkedCellService; }
  const simplemd::services::MolecularPropertiesService& getMolecularPropertiesService() {
    #if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if(_molecularPropertiesService == NULL){
      std::cout <<"ERROR coupling::interface::MDSimulation::getMolecularPropertiesService(): _molecularPropertiesService == NULL " << std::endl;
      exit(1);
    }
    #endif
    return *_molecularPropertiesService;
  }

private:
  /** @brief the coupling cell service for the coupled md simulation  */
  coupling::services::CouplingCellService<MD_DIM>* _couplingCellService;
  /** @brief bool holding the current state of the coupling: true - coupled
   * simulation and false - independent md simulation */
  bool _couplingSwitchedOn;
};
} // namespace interface
} // namespace coupling
