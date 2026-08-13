// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESPOTENTIALENERGYMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include <Kokkos_Core.hpp>

namespace simplemd {
namespace cellmappings {
class LennardJonesPotentialEnergyMapping;
}
} // namespace simplemd

/** computes the Lennard Jones energy for all atomp pairs.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::LennardJonesPotentialEnergyMapping {
public:
  LennardJonesPotentialEnergyMapping(const simplemd::services::MolecularPropertiesService& molecularPropertiesService);
  ~LennardJonesPotentialEnergyMapping() {}

  void beginCellIteration();

  void endCellIteration() {}
  KOKKOS_FUNCTION void handleCell(LinkedCell& cell) const;
  KOKKOS_FUNCTION void handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) const;

  static const bool IsParallel = true;

private:
  /** epsilon */
  const double _epsilon;
  /** sigma^6 */
  const double _sigma6;
  /** cutOffRadius*cutOffRadius */
  const double _cutOffRadiusSquared;

  /** energy to be subtracted for truncated shifted LJ potentials */
  const double _cutOffEnergy;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESPOTENTIALENERGYMAPPING_H_
