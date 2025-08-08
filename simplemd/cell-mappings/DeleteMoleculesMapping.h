// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_

#include "simplemd/LinkedCell.h"
#include <Kokkos_Core.hpp>

namespace simplemd {
namespace cellmappings {
class DeleteMoleculesMapping;
}
} // namespace simplemd

/** deletes all molecules within the cells. The molecules are not only removed from the lists in the cells.
 *  This mapping is applied in the BoundaryTreatment class: After the force computations, the ghost cells do not
 *  need to contain the molecules anymore. So, this mapping is used to delete all molecules within the ghost layers.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::DeleteMoleculesMapping {
public:
  DeleteMoleculesMapping() {}
  KOKKOS_FUNCTION ~DeleteMoleculesMapping() {}

  void beginCellIteration() {}
  void endCellIteration() {}
  KOKKOS_FUNCTION void handleCell(LinkedCell& cell) const { cell.clear(); }

  static const bool IsParallel = true;

private:
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_DELETEMOLECULESMAPPING_H_
