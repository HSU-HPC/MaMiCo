// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_

#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MoleculeService.h"
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
  ParallelBoundaryEmptyCellsMapping(simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::services::MoleculeService& moleculeService,
                                    const simplemd::services::LinkedCellService& linkedCellService)
      : _parallelTopologyService(parallelTopologyService), _moleculeService(moleculeService), _linkedCellService(linkedCellService) {}
  ~ParallelBoundaryEmptyCellsMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    // send molecules from this cell first...
    if (_parallelTopologyService.reduceGhostCellViaBuffer(cell, cellIndex, _linkedCellService, _moleculeService)) {
      // ... and erase them afterwards
      for (auto it = cell.begin(_moleculeService); it != cell.end(); it++) {
        _moleculeService.deleteMolecule(*(*it));
      }
      cell.clear(_moleculeService);
    }
  }
  static const bool IsParallel = false;

private:
  simplemd::services::ParallelTopologyService& _parallelTopologyService;
  simplemd::services::MoleculeService& _moleculeService;
  const simplemd::services::LinkedCellService& _linkedCellService;
};
#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_PARALLELBOUNDARYEMPTYCELLSMAPPING_H_
