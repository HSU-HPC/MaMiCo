// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_
#define _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_

#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/CouplingMDDefinitions.h"

namespace coupling {
  namespace paralleltopology {
    template<unsigned int dim>
    class XYZTopology;
  }
}


/** the XYZTopology orders the ranks in x-y-z manner, i.e. we obtain the rank from process coordinates (x,y,z) by
 *  z*nx*ny + y*nx + x=x + nx*(y+ny*z), where nx,ny,nz are the numbers of processes in x,y,z-direction.
 *  topologyOffset can be used to shift the whole topology by an offset of ranks.
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::paralleltopology::XYZTopology: public coupling::paralleltopology::ParallelTopology<dim> {
  public:
    XYZTopology(tarch::la::Vector<dim,unsigned int> numberProcesses,unsigned int topologyOffset):
    coupling::paralleltopology::ParallelTopology<dim>(),
    _numberProcesses(numberProcesses),
    _divisionFactor4NumberProcesses(coupling::initDivisionFactor<dim>(numberProcesses)),
    _topologyOffset(topologyOffset){}

    virtual ~XYZTopology(){}

    tarch::la::Vector<dim,unsigned int> getProcessCoordinates(unsigned int rank) const {
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      unsigned int intNumberProcesses = _numberProcesses[0]; for (unsigned int d = 1 ; d < dim; d++){ intNumberProcesses = intNumberProcesses*_numberProcesses[d]; }
      if ( (rank<_topologyOffset) || (rank>_topologyOffset+intNumberProcesses-1)){
        //TODO thrown on macroOnly services
        std::cout << "Warning coupling::paralleltopology::XYZTopology::getProcessCoordinates(): rank out of range!" << std::endl;
        std::cout << "Offset=" << _topologyOffset << ", rank=" << rank << std::endl;
      }
      std::cout << "Rank=" << rank << " corresponds to process coordinates=" << coupling::getVectorCellIndex<dim>(rank-_topologyOffset,_divisionFactor4NumberProcesses) << std::endl;
      #endif
      return coupling::getVectorCellIndex<dim>(rank-_topologyOffset,_divisionFactor4NumberProcesses);
    }


    unsigned int getRank(tarch::la::Vector<dim,unsigned int> processCoordinates) const {
      unsigned int index = processCoordinates[dim-1];
      for (int d = dim-2; d >-1; d--){
        index = _numberProcesses[d]*index + processCoordinates[d];
      }
      #if (COUPLING_MD_DEBUG==COUPLING_MD_YES)
      std::cout << "Process coordinates=" << processCoordinates << " correspond to rank=" << index+_topologyOffset << std::endl;
      #endif
      return index+_topologyOffset;
    }

  private:
    /* number of processes */
    const tarch::la::Vector<dim,unsigned int> _numberProcesses;
    /* division factor for number of processes */
    const tarch::la::Vector<dim,unsigned int> _divisionFactor4NumberProcesses;
    /* offset in ranks for linear shift of xyz-topology */
    const unsigned int _topologyOffset;
};

#endif // _MOLECULARDYNAMICS_COUPLING_PARALLELTOPOLOGY_XYZTOPOLOGY_H_

