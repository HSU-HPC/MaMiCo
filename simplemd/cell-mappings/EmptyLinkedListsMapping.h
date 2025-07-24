// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_EMPTYLINKEDLISTSMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_EMPTYLINKEDLISTSMAPPING_H_

#include "simplemd/LinkedCell.h"

namespace simplemd {
namespace cellmappings {
class EmptyLinkedListsMapping;
}
} // namespace simplemd

/** deletes all molecules from the list within the cells. The molecules are only removed from the lists in the cells,
 *  but are further available via the MoleculeService.
 *  After the force evaluation over all molecules, the DeleteMoleculesMapping is used to delete all molecules from the
 *  ghost cells. The EmptyLinkedListsMapping is used to clear the lists in the inner part of the simulation domain.
 *  Here, we still need to store the molecules but need to sort these molecules again into the correct cells, as soon
 *  as the time integration is carried out.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::EmptyLinkedListsMapping {
public:
  EmptyLinkedListsMapping(simplemd::services::MoleculeService& moleculeService) : _moleculeService(moleculeService) {}
  ~EmptyLinkedListsMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}
  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) { cell.clear(_moleculeService); }

private:
  simplemd::services::MoleculeService& _moleculeService;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_EMPTYLINKEDLISTSMAPPING_H_
