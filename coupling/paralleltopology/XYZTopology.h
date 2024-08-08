// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/paralleltopology/ParallelTopology.h"

namespace coupling {
namespace paralleltopology {
template <unsigned int dim> class XYZTopology;
}
} // namespace coupling

/** the XYZTopology orders the ranks in x-y-z manner, i.e. we obtain the rank
 * from process coordinates (x,y,z) by z*nx*ny + y*nx + x=x + nx*(y+ny*z), where
 * nx,ny,nz are the numbers of processes in x,y,z-direction. topologyOffset can
 * be used to shift the whole topology by an offset of ranks. Derived class from
 * the class ParallelTopology. E.g. assuming ParallelTopologyType = XYZ and there
 * is a cubic domain, splitted into 8 sub-domains (2 sub-domains in each
 * dimension). Then the ordering of the MPI processes is: Rank=0 for x=0,y=0,z=0.
 * Rank=1 for x=1,y=0,z=0. Rank=2 for x=0,y=1,z=0. Rank=3 for x=1,y=1,z=0. Rank=4
 * for x=0,y=0,z=1. Rank=5 for x=1,y=0,z=1. Rank=6 for x=0,y=1,z=1. Rank=7 for
 * x=1,y=1,z=1.
 * @brief The XYZTopology orders the ranks in x-y-z manner.
 * @tparam dim Number of dimensions; it can be 1, 2 or 3
 * @author Philipp Neumann
 * @todo Philipp could you please take a look on this class
 */
template <unsigned int dim> class coupling::paralleltopology::XYZTopology : public coupling::paralleltopology::ParallelTopology<dim> {
public:
  /** Constructor */
  XYZTopology(tarch::la::Vector<dim, unsigned int> numberProcesses)
      : coupling::paralleltopology::ParallelTopology<dim>(), _numberProcesses(numberProcesses),
        _divisionFactor4NumberProcesses(coupling::initDivisionFactor<dim>(numberProcesses)) {}

  /** Destructor */
  virtual ~XYZTopology() {}

  tarch::la::Vector<dim, unsigned int> getProcessCoordinates(unsigned int rank, unsigned int topologyOffset) const {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Rank=" << rank
              << " corresponds to process coordinates=" << coupling::getVectorCellIndex<dim>(rank - topologyOffset, _divisionFactor4NumberProcesses)
              << std::endl;
#endif
    return coupling::getVectorCellIndex<dim>(rank - topologyOffset, _divisionFactor4NumberProcesses);
  }

  unsigned int getRank(tarch::la::Vector<dim, unsigned int> processCoordinates, unsigned int topologyOffset) const {
    unsigned int index = processCoordinates[dim - 1];
    for (int d = dim - 2; d > -1; d--) {
      index = _numberProcesses[d] * index + processCoordinates[d];
    }
    return index + topologyOffset;
  }

private:
  /* number of processes */
  const tarch::la::Vector<dim, unsigned int> _numberProcesses;
  /* division factor for number of processes */
  const tarch::la::Vector<dim, unsigned int> _divisionFactor4NumberProcesses;
};

#endif // _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_
