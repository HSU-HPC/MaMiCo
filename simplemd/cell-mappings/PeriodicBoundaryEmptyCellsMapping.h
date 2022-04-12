// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICBOUNDARYEMPTYCELLSMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICBOUNDARYEMPTYCELLSMAPPING_H_

#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/services/ParallelTopologyService.h"

namespace simplemd {
namespace cellmappings {
class PeriodicBoundaryEmptyCellsMapping;
}
} // namespace simplemd

/** this class is used to handle periodic boundaries. It takes all particles in
 * outer (ghost) cells, puts them into the respective inner cells (w.r.t.
 * periodicity) and deletes them from the ghost cell. In case of (MPI-)parallel
 * computations, the molecules are automatically adapted in the positions and
 * sent to the respective neighboring process (w.r.t. periodicity).
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping {
public:
  PeriodicBoundaryEmptyCellsMapping(simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::services::MoleculeService& moleculeService,
                                    simplemd::services::LinkedCellService& linkedCellService);
  ~PeriodicBoundaryEmptyCellsMapping() {}

  /** sets the global domain size (hopefully received from the
   * ParallelTopologyService...) */
  void setDomainSize(const tarch::la::Vector<MD_DIM, double>& domainSize);

  void setProcessCoordinates(const tarch::la::Vector<MD_DIM, unsigned int>& processCoordinates);

  void setNumberOfProcesses(const tarch::la::Vector<MD_DIM, unsigned int>& numberProcesses);

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell& cell, const unsigned int& cellIndex);

private:
  simplemd::services::ParallelTopologyService& _parallelTopologyService;
  simplemd::services::MoleculeService& _moleculeService;
  simplemd::services::LinkedCellService& _linkedCellService;
  /** domain size and offset */
  tarch::la::Vector<MD_DIM, double> _domainSize;
  tarch::la::Vector<MD_DIM, unsigned int> _processCoordinates;
  tarch::la::Vector<MD_DIM, unsigned int> _numberProcesses;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_PERIODICBOUNDARYEMPTYCELLSMAPPING_H_
