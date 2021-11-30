// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_

#include "coupling/ParticleInsertion.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class NoParticleInsertion;
}

/** @brief empty implementation of particle insertion (no particle insertion/removal).
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2, or 3
 */
template<class LinkedCell,unsigned int dim>
class coupling::NoParticleInsertion: public coupling::ParticleInsertion<LinkedCell,dim> {
  public:
    /** @brief dummy function, doesn't change or do something
     *  @returns coupling::ParticleInsertion::NoAction since no mass is deleted or inserted */
    virtual typename coupling::ParticleInsertion<LinkedCell,dim>::Action insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      const tarch::la::Vector<dim,double>& meanVelocity,
      const double &temperature,
      const coupling::BoundaryForceController<LinkedCell,dim>& boundaryForceController
    ) {return coupling::ParticleInsertion<LinkedCell,dim>::NoAction;}
    
    /** @brief returns false, since no energy landscape is needed in the dummy class
    *  @returns a bool that indicates if a potential energy landscape is neccessary (true) or not (false) */
    virtual bool requiresPotentialEnergyLandscape(){return false;}

    /** @brief a simple constructor*/
    NoParticleInsertion(): coupling::ParticleInsertion<LinkedCell,dim>(1){}

    /** @brief a simple destructor*/
    virtual ~NoParticleInsertion(){}

};
#endif // _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_
