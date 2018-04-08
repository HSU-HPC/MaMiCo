// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MolecularPropertiesService.h"
#include "simplemd/services/ExternalForceService.h"
#include "tarch/la/ScalarOperations.h"

namespace simplemd {
  namespace cellmappings {
    class LennardJonesForceMapping;
  }
}


/** applies the Lennard-Jones force to all particle pairs.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::LennardJonesForceMapping {
  public:
    LennardJonesForceMapping(simplemd::services::ExternalForceService &externalForceService,const simplemd::services::MolecularPropertiesService &molecularPropertiesService);
    ~LennardJonesForceMapping(){}

    void beginCellIteration();

    void endCellIteration(){}
    void handleCell(const LinkedCell& cell,const unsigned int& cellIndex) const;
    void handleCellPair(const LinkedCell& cell1, const LinkedCell& cell2,const unsigned int& cellIndex1, const unsigned int& cellIndex2) const;

    /** returns the force acting on a particle placed at position1, resulting from an interaction of the particles at
     *  positions position1 and position2. Remark: The force on the particle at position2 is just (-1.0)*returnValue.
     */
    tarch::la::Vector<MD_DIM,double> getLennardJonesForce(
      const tarch::la::Vector<MD_DIM,double>& position1,const tarch::la::Vector<MD_DIM,double>& position2
    ) const;

  private:
    /** epsilon */
    const double _epsilon;
    /** sigma^6 */
    const double _sigma6;
    /** cutOffRadius*cutOffRadius */
    const double _cutOffRadiusSquared;
    /** external forces*/
    simplemd::services::ExternalForceService &_externalForceService;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_LENNARDJONESFORCEMAPPING_H_

