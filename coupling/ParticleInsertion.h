// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_

#include <iostream>
#include <cstdlib>
#include "coupling/BoundaryForceController.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class ParticleInsertion;
}


/** interface for particle insertion/removal on macroscopic cell basis.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::ParticleInsertion {
  public:
    /** this state is returned by the insertDeleteMass() function and tells the user, if mass was inserted/deleted
     *  or if nothing happened at all.
     */
    enum Action{NoAction=0,Insertion=1,Deletion=2};

    virtual typename coupling::ParticleInsertion<LinkedCell,dim>::Action insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      const tarch::la::Vector<dim,double>& meanVelocity,
      const double &temperature,
      bool toBePlotted,
      const coupling::BoundaryForceController<LinkedCell,dim>& boundaryForceController
    )  = 0;

    /** returns true, if the particle insertion requires information on the potential energy landscape. The USHER
     *  algorithm requires it. Other algorithms may not require it; the trivial NoParticleInsertion-implementation
     *  which does not do anything obviously does not require a valid potential energy landscape.
     */
    virtual bool requiresPotentialEnergyLandscape() = 0;

    ParticleInsertion(unsigned int insertDeleteMassEveryTimestep): _insertDeleteMassEveryTimestep(insertDeleteMassEveryTimestep){
      if (_insertDeleteMassEveryTimestep==0){
        std::cout << "ERROR ParticleInsertion::ParticleInsertion(..): _insertDeleteMassEveryTimestep=0!" << std::endl; exit(EXIT_FAILURE);
      }
    }
    virtual ~ParticleInsertion(){}

    bool insertDeleteMassAtTimestep(unsigned int t) const {
      return (t%_insertDeleteMassEveryTimestep==0);
    }

  private:
    const unsigned int _insertDeleteMassEveryTimestep;
};
#endif // _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
