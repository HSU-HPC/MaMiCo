// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/ExternalForceService.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "tarch/la/ScalarOperations.h"

#include <Kokkos_Core.hpp>

namespace simplemd {
namespace cellmappings {
class LennardJonesForceMapping;
}
} // namespace simplemd

/** applies the Lennard-Jones force to all particle pairs.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::LennardJonesForceMapping {
public:
  LennardJonesForceMapping(simplemd::services::ExternalForceService& externalForceService,
                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService);
  KOKKOS_FUNCTION ~LennardJonesForceMapping() {}

  KOKKOS_FUNCTION void beginCellIteration();

  KOKKOS_FUNCTION void endCellIteration() {}
  KOKKOS_FUNCTION void handleCell(LinkedCell& cell) const;
  KOKKOS_FUNCTION void handleCellPair(const LinkedCell& cell1, const LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) const;

  /** returns the force acting on a particle placed at position1, resulting from an interaction of the particles at
   *  positions position1 and position2. Remark: The force on the particle at position2 is just (-1.0)*returnValue.
   */
  KOKKOS_FUNCTION tarch::la::Vector<MD_DIM, double> getLennardJonesForce(const tarch::la::Vector<MD_DIM, double>& position1,
                                                         const tarch::la::Vector<MD_DIM, double>& position2) const;

  static const bool IsParallel = true;

private:
  /** epsilon */
  const double _epsilon;
  /** sigma^6 */
  const double _sigma6;
  /** cutOffRadius*cutOffRadius */
  const double _cutOffRadiusSquared;
  /** external forces*/
  tarch::la::Vector<MD_DIM, double> _externalForce;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_
