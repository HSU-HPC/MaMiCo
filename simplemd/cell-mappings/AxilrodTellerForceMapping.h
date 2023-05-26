// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERFORCEMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERFORCEMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/ExternalForceService.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "tarch/la/ScalarOperations.h"
#include <vector>

namespace simplemd {
namespace cellmappings {
class AxilrodTellerForceMapping;
}
} // namespace simplemd

/** applies the Axilrod-Teller force to all particle triplets.
 *  @author Maximilian Mayr
 */
class simplemd::cellmappings::AxilrodTellerForceMapping {
public:
  AxilrodTellerForceMapping(simplemd::services::ExternalForceService& externalForceService,
                           const simplemd::services::MolecularPropertiesService& molecularPropertiesService);
  ~AxilrodTellerForceMapping() {}

  void beginCellIteration();

  void endCellIteration() {}
  void handleCell(const LinkedCell& cell, const unsigned int& cellIndex) const;
  void handleCellPair(const LinkedCell& cell1, const LinkedCell& cell2, const unsigned int& cellIndex1, const unsigned int& cellIndex2) const;
  void handleCellTriplet(const LinkedCell& cell1, const LinkedCell& cell2, const LinkedCell& cell3,
                         const unsigned int& cellIndex1, const unsigned int& cellIndex2, const unsigned int& cellIndex3) const;

  /** returns the force acting on the particles placed at position1, position2 and position3 resulting from an interaction of the particles at
   *  positions position1, position2 and position3. 
   */
  std::vector<tarch::la::Vector<MD_DIM, double>> getAxilrodTellerForce(const tarch::la::Vector<MD_DIM, double>& position1,
                                                         const tarch::la::Vector<MD_DIM, double>& position2,
                                                         const tarch::la::Vector<MD_DIM, double>& position3) const;

private:
  /** v */
  const double _v;
  /** cutOffRadius*cutOffRadius */
  const double _cutOffRadiusSquared;
  /** external forces*/
  simplemd::services::ExternalForceService& _externalForceService;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_AXILRODTELLERFORCEMAPPING_H_