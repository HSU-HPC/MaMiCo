// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_

#include "simplemd/MolecularDynamicsSimulation.h"
#include "coupling/services/MacroscopicCellService.h"

namespace coupling {
namespace solvers { class CoupledMolecularDynamicsSimulation; }
}

/** simulation class for coupled MD simulations; thus, the implementation of one
 * timestep slightly
 *  differs from the one of the base class.
 *
 *  @author Philipp Neumann
 */
class coupling::solvers::CoupledMolecularDynamicsSimulation
    : public simplemd::MolecularDynamicsSimulation {
public:
  CoupledMolecularDynamicsSimulation(
      const simplemd::configurations::MolecularDynamicsConfiguration &
          configuration);
  virtual ~CoupledMolecularDynamicsSimulation() {}

  void simulateOneCouplingTimestep(const unsigned int &t);

  void switchOnCoupling() { _couplingSwitchedOn = true; }
  void switchOffCoupling() { _couplingSwitchedOn = false; }

  void setMacroscopicCellService(coupling::services::MacroscopicCellService<
      MD_DIM> *macroscopicCellService) {
    _macroscopicCellService = macroscopicCellService;
  }

  /** returns the boundary treatment of the simulation; this is needed by the
   * coupling to synchronize molecules
   *  in boundary regions and on different processes after mass insertion/
   * deletion.
   */
  simplemd::BoundaryTreatment &getBoundaryTreatment() {
    return *_boundaryTreatment;
  }

  /** returns the parallel topology service; this is needed by the MD solver
   * interface in coupling. */
  simplemd::services::ParallelTopologyService &getParallelTopologyService() {
    return *_parallelTopologyService;
  }
  /** returns the molecule service; this is needed by the MD solver interface in
   * coupling. */
  simplemd::services::MoleculeService &getMoleculeService() {
    return *_moleculeService;
  }
  simplemd::services::LinkedCellService &getLinkedCellService() {
    return *_linkedCellService;
  }
  const simplemd::services::MolecularPropertiesService &
  getMolecularPropertiesService() {
    return *_molecularPropertiesService;
  }

private:
  coupling::services::MacroscopicCellService<MD_DIM> *_macroscopicCellService;
  bool _couplingSwitchedOn;
};
#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUPLEDMOLECULARDYNAMICSSIMULATION_H_
