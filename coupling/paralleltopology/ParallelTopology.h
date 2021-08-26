// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_

#include "tarch/la/Vector.h"

namespace coupling {
  namespace paralleltopology {
    template<unsigned int dim>
    class ParallelTopology;
 }
}


/** interface for different parallel topologies. This can be adapted to the respective Cartesian topology
 *  that is applied to the MD simulation. This class solely performs the conversion rank<-> process coordinates,
 *  assuming a Cartesian grid-like splitting of the MD domain.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::paralleltopology::ParallelTopology {
  public:
    virtual ~ParallelTopology(){}
    /** converts process coordinates into a rank */
    virtual unsigned int getRank(tarch::la::Vector<dim,unsigned int> processCoordinates) const = 0;
    /** converts rank into process coordinates */
    virtual tarch::la::Vector<dim,unsigned int> getProcessCoordinates(unsigned int rank) const = 0;
};

#endif // _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_
