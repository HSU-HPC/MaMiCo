// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_

#include "tarch/la/Vector.h"

namespace coupling {
	/** namespace coupling */
  namespace paralleltopology {
	  /** namespace paralleltopology */
    template<unsigned int dim>
    class ParallelTopology;
 }
}


/** interface for different parallel topologies. This can be adapted to the respective Cartesian topology
 *  that is applied to the MD simulation. This class solely performs the conversion rank<-> process coordinates,
 *  assuming a Cartesian grid-like splitting of the MD domain. It is the base class of the classe ZYXTopology and XYZTopology.
 *	@brief This class performs the conversion rank <-> process coordinates.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::paralleltopology::ParallelTopology {
  public:
	/** Destructor */
    virtual ~ParallelTopology(){}
    /** This function takes process coordinates as input and returns the correpsponding rank.
	 *	@brief converts process coordinates into a rank.
	 *	@param processCoordinates Process coordinates
	 */
    virtual unsigned int getRank(tarch::la::Vector<dim,unsigned int> processCoordinates) const = 0;
    /** This function takes rank as input and return the correpsponding process coordinates.
	 *	@brief converts rank into process coordinates.
	 *	@param rank Rank
	 */
    virtual tarch::la::Vector<dim,unsigned int> getProcessCoordinates(unsigned int rank) const = 0;
};

#endif // _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGY_H_
