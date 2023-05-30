// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERPOTENTIALENERGYMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "tarch/la/ScalarOperations.h"

namespace simplemd {
namespace cellmappings {
class AxilrodTellerPotentialEnergyMapping;
}
} // namespace simplemd

/** computes the Axilrod-Teller energy for all atomp triplets.
 *  @author Maximilian Mayr
 */
class simplemd::cellmappings::AxilrodTellerPotentialEnergyMapping {
public:
  AxilrodTellerPotentialEnergyMapping(const simplemd::services::MolecularPropertiesService& molecularPropertiesService);
  ~AxilrodTellerPotentialEnergyMapping() {}

  void beginCellIteration();

  void endCellIteration() {}
  void handleCell(LinkedCell& cell, const unsigned int& cellIndex);
  void handleCellPair(LinkedCell& cell1, LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2);
  void handleCellTriplet(LinkedCell& cell1, LinkedCell& cell2, LinkedCell& cell3,
                         const unsigned int& cellIndex1, const unsigned int& cellIndex2, const unsigned int& cellIndex3);
  tarch::la::Vector<3, double> getAxilrodTellerPotentialEnergy(const tarch::la::Vector<MD_DIM, double>& position1,
                        const tarch::la::Vector<MD_DIM, double>& position2, const tarch::la::Vector<MD_DIM, double>& position3);

private:
  /** v */
  const double _v;
  /** cutOffRadius*cutOffRadius */
  const double _cutOffRadiusSquared;

  /** energy to be subtracted for truncated shifted LJ potentials */
  const double _cutOffEnergy;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERPOTENTIALENERGYMAPPING_H_