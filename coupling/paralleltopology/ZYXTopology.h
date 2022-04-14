// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_ZYXTOPOLOGY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_ZYXTOPOLOGY_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/paralleltopology/ParallelTopology.h"

namespace coupling {
namespace paralleltopology {
template <unsigned int dim> class ZYXTopology;
}
} // namespace coupling

/** In the ZYXTopology, the process coordinates convert to a rank as
 *  rank = x*ny*nz + y*nz + z = z + nz*(y+ny*x) (for 3D).
 *  topologyOffset is used for linearized access of multiple MD instances.
 *Derived class from the class ParallelTopology. E.g. assuming
 *ParallelTopologyType = XYZ and there is a cubic domain, splitted into 8
 *sub-domains (2 sub-domains in each dimension). Then the ordring of the MPI
 *processes is: Rank=0 for x=0,y=0,z=0. Rank=1 for x=1,y=0,z=0. Rank=2 for
 *x=0,y=1,z=0. Rank=3 for x=1,y=1,z=0. Rank=4 for x=0,y=0,z=1. Rank=5 for
 *x=1,y=0,z=1. Rank=6 for x=0,y=1,z=1. Rank=7 for x=1,y=1,z=1.
 *	@brief The XYZTopology orders the ranks in z-y-x manner.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 *	@todo Philipp could you please take a look on this class
 */
template <unsigned int dim> class coupling::paralleltopology::ZYXTopology : public coupling::paralleltopology::ParallelTopology<dim> {
public:
  /** Constructor */
  ZYXTopology(tarch::la::Vector<dim, unsigned int> numberProcesses, unsigned int topologyOffset)
      : coupling::paralleltopology::ParallelTopology<dim>(), _numberProcesses(numberProcesses),
        _divisionFactor4NumberProcesses(initDivisionFactor(numberProcesses)), _topologyOffset(topologyOffset) {}

  /** Destructor */
  virtual ~ZYXTopology() {}

  tarch::la::Vector<dim, unsigned int> getProcessCoordinates(unsigned int rank) const {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    unsigned int intNumberProcesses = _numberProcesses[0];
    for (unsigned int d = 1; d < dim; d++) {
      intNumberProcesses = intNumberProcesses * _numberProcesses[d];
    }
    if ((rank < _topologyOffset) || (rank > _topologyOffset + intNumberProcesses - 1)) {
      std::cout << "Warning "
                   "coupling::paralleltopology::ZYXTopology::"
                   "getProcessCoordinates(): rank out of range!"
                << std::endl;
      std::cout << "Offset=" << _topologyOffset << ", rank=" << rank << std::endl;
    }
#endif
    tarch::la::Vector<dim, unsigned int> processCoordinates(0);
    unsigned int help = rank - _topologyOffset;
    for (unsigned int d = 0; d < dim; d++) {
      processCoordinates[d] = help / _divisionFactor4NumberProcesses[d];
      help = help - processCoordinates[d] * _divisionFactor4NumberProcesses[d];
    }
    return processCoordinates;
  }

  /** computes the rank as shown above, see second formula of class definition
   */
  unsigned int getRank(tarch::la::Vector<dim, unsigned int> processCoordinates) const {
    unsigned int rank = processCoordinates[0];
    for (unsigned int d = 1; d < dim; d++) {
      rank = rank * _numberProcesses[d] + processCoordinates[d];
    }
    return rank + _topologyOffset;
  }

private:
  /** sets the division factor for each vector entry. For ZYX, this corresponds
   * to (in 3D) (ny*nz,nz,1) and to (2D) (ny,1). */
  tarch::la::Vector<dim, unsigned int> initDivisionFactor(tarch::la::Vector<dim, unsigned int> numberProcesses) const {
    tarch::la::Vector<dim, unsigned int> div(1);
    for (int d = dim - 2; d > -1; d--) {
      div[d] = div[d + 1] * numberProcesses[d + 1];
    }
    return div;
  }

  /* number of processes */
  const tarch::la::Vector<dim, unsigned int> _numberProcesses;
  /* division factor */
  const tarch::la::Vector<dim, unsigned int> _divisionFactor4NumberProcesses;
  const unsigned int _topologyOffset;
};
#endif // _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_ZYXTOPOLOGY_H_
