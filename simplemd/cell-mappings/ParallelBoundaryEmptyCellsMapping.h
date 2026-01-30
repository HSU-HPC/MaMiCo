// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_

#include "simplemd/services/ParallelTopologyService.h"

namespace simplemd {
namespace cellmappings {
class ParallelBoundaryEmptyCellsMapping;
}
} // namespace simplemd

/** shall only be applied to PARALLEL_BOUNDARY ghost cells. All particles of each cell are sent to the
 *  respective neighbour with the inner cell corresponding to these ghost boundary cells.
 *  Afterwards, the list of the ghost cell is cleared and the molecules are deleted from the ghost cell.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::ParallelBoundaryEmptyCellsMapping {
public:
  ParallelBoundaryEmptyCellsMapping(simplemd::services::ParallelTopologyService& parallelTopologyService, const simplemd::MoleculeContainer& moleculeContainer)
      : _parallelTopologyService(parallelTopologyService), _moleculeContainer(moleculeContainer) {}
  ~ParallelBoundaryEmptyCellsMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell& cell) {
    // send molecules from this cell first...
    if (_parallelTopologyService.reduceGhostCellViaBuffer(cell, cell.getIndex(), _moleculeContainer)) {
      // ... and erase them afterwards
      cell.clear();
    }
  }
  static const bool IsParallel = false;

private:
  simplemd::services::ParallelTopologyService& _parallelTopologyService;
  const simplemd::MoleculeContainer& _moleculeContainer;
};
#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_
