// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULEMAPPINGS_UPDATELINKEDCELLLISTSMAPPING_H_
#define _MOLECULARDYNAMICS_MOLECULEMAPPINGS_UPDATELINKEDCELLLISTSMAPPING_H_

#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/ParallelTopologyService.h"

namespace simplemd {
namespace moleculemappings {
class UpdateLinkedCellListsMapping;
}
} // namespace simplemd

/** cell mapping used to sort molecules into the respective linked cells.
 *  This step is done twice in the algorithm:
 *  - Once at the startup of the simulation. When the LinkedCellService is
 * initialised, the molecules need to be sorted into the corresponding cells.
 *  - After each time integration step for the molecules, the cells are emptied
 * and refilled by this mapping.
 *
 *  @author Philipp Neumann
 */
class simplemd::moleculemappings::UpdateLinkedCellListsMapping {
public:
  UpdateLinkedCellListsMapping(const simplemd::services::ParallelTopologyService& parallelTopologyService,
                               simplemd::services::LinkedCellService& linkedCellService, simplemd::services::MoleculeService& moleculeService)
      : _parallelTopologyService(parallelTopologyService), _linkedCellService(linkedCellService), _moleculeService(moleculeService) {}
  ~UpdateLinkedCellListsMapping() {}

  void beginMoleculeIteration();
  void endMoleculeIteration() {}
  void handleMolecule(Molecule& molecule);

private:
  const simplemd::services::ParallelTopologyService& _parallelTopologyService;
  simplemd::services::LinkedCellService& _linkedCellService;
  simplemd::services::MoleculeService& _moleculeService;
#if (MD_ERROR == MD_YES)
  /** global domain size */
  tarch::la::Vector<MD_DIM, double> _domainSize;
#endif

  /** global domain offset */
  tarch::la::Vector<MD_DIM, double> _domainOffset;

  /** mesh width */
  tarch::la::Vector<MD_DIM, double> _meshWidth;

  /** global index of the first cell of this domain */
  tarch::la::Vector<MD_DIM, int> _globalIndexOfFirstCell;

  /** local index of the first cell within this domain */
  tarch::la::Vector<MD_DIM, int> _localIndexOfFirstCell;
};

#endif // _MOLECULARDYNAMICS_MOLECULEMAPPINGS_UPDATELINKEDCELLLISTSMAPPING_H_
