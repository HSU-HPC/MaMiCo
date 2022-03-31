// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICANDPARALLELBOUNDARYFILLCELLSMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICANDPARALLELBOUNDARYFILLCELLSMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/MolecularDynamicsDefinitions.h"
#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/services/ParallelTopologyService.h"

namespace simplemd {
namespace cellmappings {
class PeriodicAndParallelBoundaryFillCellsMapping;
}
} // namespace simplemd

/** This mapping is used at the very beginning of each MD timestep. In order to
 * populate the ghost cells of
 *  periodic and parallel boundaries, this mapping triggers a broadcasting of
 * inner cell information to the respective
 *  parallel and periodic ghost cells. At this, the molecules are either
 * directly manipulated and sent to the
 *  respective process, or their position is locally adapted according to the
 * periodic boundary condition, and they
 *  are cloned and put into the corresponding local ghost cells.
 *  This mapping is triggered from the BoundaryTreatment class.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::PeriodicAndParallelBoundaryFillCellsMapping {
public:
  PeriodicAndParallelBoundaryFillCellsMapping(simplemd::services::ParallelTopologyService &parallelTopologyService,
                                              simplemd::services::MoleculeService &moleculeService, simplemd::services::LinkedCellService &linkedCellService)
      : _parallelTopologyService(parallelTopologyService), _moleculeService(moleculeService), _linkedCellService(linkedCellService) {}
  ~PeriodicAndParallelBoundaryFillCellsMapping() {}

  void setDomainSize(const tarch::la::Vector<MD_DIM, double> &domainSize) { _domainSize = domainSize; }

  void beginCellIteration() {}
  void endCellIteration() {}
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex);

private:
  simplemd::services::ParallelTopologyService &_parallelTopologyService;
  simplemd::services::MoleculeService &_moleculeService;
  simplemd::services::LinkedCellService &_linkedCellService;
  /** domain size */
  tarch::la::Vector<MD_DIM, double> _domainSize;
};
#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICANDPARALLELBOUNDARYFILLCELLSMAPPING_H_
