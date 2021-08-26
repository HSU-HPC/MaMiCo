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
    /** parallel topology types that are supported */
    enum ParallelTopologyType{UNDEFINED=0,XYZ=1,ZYX=2};
  }
}


/** creates the parallel topology from a given topology type and a number of processes (typically read from a configuration).
 *  @author Philipp Neumann
 */
class coupling::paralleltopology::ParallelTopologyFactory {
  public:
    template<unsigned int dim>
    static coupling::paralleltopology::ParallelTopology<dim>* getParallelTopology(
      coupling::paralleltopology::ParallelTopologyType type, tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int topologyOffset
    ) {
      if (type == XYZ){
        return new coupling::paralleltopology::XYZTopology<dim>(numberProcesses,topologyOffset);
      } else if (type == ZYX){
        return new coupling::paralleltopology::ZYXTopology<dim>(numberProcesses,topologyOffset);
      } else {
        return NULL;
      }
    }
};
#endif
