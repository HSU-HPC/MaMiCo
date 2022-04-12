// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_

#include "coupling/services/MacroscopicCellService.h"
#include "simplemd/MolecularDynamicsSimulation.h"

namespace coupling {
namespace solvers {
class CoupledMolecularDynamicsSimulation;
}
} // namespace coupling

/** simulation class for coupled MD simulations; thus, the implementation of one
 * timestep slightly differs from the one of the base class.
 *  @author Philipp Neumann
 *  @todo since simpleMD is not part of the docu, somthing is missing here?!?
 * How to deal with?
 */
class coupling::solvers::CoupledMolecularDynamicsSimulation : public simplemd::MolecularDynamicsSimulation {
public:
  /** @brief a simple constructor
   *  @param configuration configuration information for the md setup*/
  CoupledMolecularDynamicsSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration);

  /** @brief a dummy destructor */
  virtual ~CoupledMolecularDynamicsSimulation() {}

  /** @brief simulates one coupled time step of the md solver */
  void simulateOneCouplingTimestep(const unsigned int& t);

  /** @brief enables the coupling, e.g. sets the boundary condition accordingly
   */
  void switchOnCoupling() { _couplingSwitchedOn = true; }

  /** @brief turns off the coupling parts within md, md runs according to the
   * basic md method */
  void switchOffCoupling() { _couplingSwitchedOn = false; }

  /** @brief set a macroscopic cell service for the coupled md simulation
   *  @param macroscopicCellService the macroscopicCellService to set */
  void setMacroscopicCellService(coupling::services::MacroscopicCellService<MD_DIM>* macroscopicCellService) {
    _macroscopicCellService = macroscopicCellService;
  }

  /** this is needed by the coupling to synchronize molecules
   *  in boundary regions and on different processes after mass insertion/
   * deletion.
   *  @brief returns the boundary treatment of the simulation;
   *  @returns the boundary treatment of the simulation */
  simplemd::BoundaryTreatment& getBoundaryTreatment() { return *_boundaryTreatment; }

  /** this is needed by the MD solver interface in coupling.
   *  @brief returns the parallel topology service;
   *  @returns the parallel topology service*/
  simplemd::services::ParallelTopologyService& getParallelTopologyService() { return *_parallelTopologyService; }

  /** this is needed by the MD solver interface in coupling.
   *  @brief returns the molecule service;
   *  @returns the molecule service */
  simplemd::services::MoleculeService& getMoleculeService() { return *_moleculeService; }

  /** @brief getter for the LinkedCellService
   *  @returns the LinkedCellService */
  simplemd::services::LinkedCellService& getLinkedCellService() { return *_linkedCellService; }

  /** @brief getter for the molecular properties service
   *  @returns the molecular properties service */
  const simplemd::services::MolecularPropertiesService& getMolecularPropertiesService() { return *_molecularPropertiesService; }

private:
  /** @brief the macroscopic cell service for the coupled md simulation  */
  coupling::services::MacroscopicCellService<MD_DIM>* _macroscopicCellService;
  /** @brief bool holding the current state of the coupling: true - coupled
   * simulation and false - independent md simulation */
  bool _couplingSwitchedOn;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_
