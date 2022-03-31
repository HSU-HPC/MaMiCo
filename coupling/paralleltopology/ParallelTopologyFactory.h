// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGYFACTORY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_PARALLELTOPOLOGYFACTORY_H_

#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/XYZTopology.h"
#include "coupling/paralleltopology/ZYXTopology.h"

namespace coupling {
namespace paralleltopology {
class ParallelTopologyFactory;
/** parallel topology types that are supported.
	 *	@enum ParallelTopologyType specifies the ordering of the MPI processes (for
domain decomposition).
	 *	As an example: assuming ParallelTopologyType = XYZ and there is a cubic
domain, splitted into 8 sub-domains (2 sub-domains in each dimension). Then the
ordring of the MPI processes is:
	 *	Rank=0 for x=0,y=0,z=0. Rank=1 for x=0,y=0,z=1. Rank=2 for x=0,y=1,z=0.
Rank=3 for x=0,1=0,z=1.
	 *	Rank=4 for x=1,y=0,z=0. Rank=5 for x=1,y=0,z=1. Rank=6 for x=1,y=1,z=0.
Rank=7 for x=1,y=1,z=1.
 	 */
enum ParallelTopologyType {
  UNDEFINED = 0 /**< gherghere gher */
      ,
  XYZ = 1 /**< the XYZTopology orders the ranks in x-y-z manner, i.e. we obtain
   the rank from process coordinates (x,y,z) by z*nx*ny + y*nx + x=x +
   nx*(y+ny*z)*/
      ,
  ZYX = 2 /**< the ZYXTopology orders the ranks in z-y-x manner, i.e. we obtain
   the rank from process coordinates (z,y,x) by x*nz*ny + y*nz + z=z +
   nz*(y+ny*x)*/
};
}
}
/** This class creates the parallel topology from a given topology type and a
 * number of processes (typically read from a configuration). assuming
 * ParallelTopologyType = XYZ and there is a cubic domain, splitted into 8
 * sub-domains (2 sub-domains in each dimension). Then the ordring of the MPI
 * processes is:
 *	Rank=0 for x=0,y=0,z=0. Rank=1 for x=0,y=0,z=1. Rank=2 for x=0,y=1,z=0.
 * Rank=3 for x=0,1=0,z=1.
 *	Rank=4 for x=1,y=0,z=0. Rank=5 for x=1,y=0,z=1. Rank=6 for x=1,y=1,z=0.
 * Rank=7 for x=1,y=1,z=1.
 *	@brief creates the parallel topology from a given topology type and a number
 * of processes
 *  @author Philipp Neumann
 */
class coupling::paralleltopology::ParallelTopologyFactory {
public:
  /** @brief This template function takes ParallelTopologyType and the number of
processes as inputs and returns a pointer to the created parallel topology.
	 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
	 *	@param type ParallelTopologyType (rdering of the MPI processes)
	 *	@param numberProcesses number of processes
	 */
  template <unsigned int dim>
  static coupling::paralleltopology::ParallelTopology<dim> *
  getParallelTopology(coupling::paralleltopology::ParallelTopologyType type,
                      tarch::la::Vector<dim, unsigned int> numberProcesses,
                      unsigned int topologyOffset) {
    if (type == XYZ) {
      return new coupling::paralleltopology::XYZTopology<dim>(numberProcesses,
                                                              topologyOffset);
    } else if (type == ZYX) {
      return new coupling::paralleltopology::ZYXTopology<dim>(numberProcesses,
                                                              topologyOffset);
    } else {
      return NULL;
    }
  }
};
#endif
