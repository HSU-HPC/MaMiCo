// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MoleculeService.h"

namespace simplemd {
namespace cellmappings {
class DeleteMoleculesMapping;
}
} // namespace simplemd

/** deletes all molecules within the cells. The molecules are not only removed
 * from the lists in the cells, they are also deleted from the MoleculeService,
 * i.e. they are completely removed. This mapping is applied in the
 * BoundaryTreatment class: After the force computations, the ghost cells do not
 *  need to contain the molecules anymore. So, this mapping is used to delete
 * all molecules within the ghost layers.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::DeleteMoleculesMapping {
public:
  DeleteMoleculesMapping(simplemd::services::MoleculeService &moleculeService)
      : _moleculeService(moleculeService) {}
  ~DeleteMoleculesMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex);

private:
  simplemd::services::MoleculeService &_moleculeService;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_
