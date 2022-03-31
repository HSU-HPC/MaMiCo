// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_

#include "coupling/ParticleInsertion.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class NoParticleInsertion;
}

/** empty implementation of particle insertion (no particle insertion/removal).
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::NoParticleInsertion
    : public coupling::ParticleInsertion<LinkedCell, dim> {
public:

  virtual typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,
                                                               dim> &cell,
      const tarch::la::Vector<dim, double> &macroscopicCellPosition,
      const tarch::la::Vector<dim, double> &macroscopicCellSize,
      const tarch::la::Vector<dim, double> &meanVelocity,
      const double &temperature,
      const coupling::BoundaryForceController<LinkedCell, dim> &
          boundaryForceController) {
    return coupling::ParticleInsertion<LinkedCell, dim>::NoAction;
  }

  virtual bool requiresPotentialEnergyLandscape() { return false; }

  NoParticleInsertion() : coupling::ParticleInsertion<LinkedCell, dim>(1) {}
  virtual ~NoParticleInsertion() {}

};
#endif // _MOLECULARDYNAMICS_COUPLING_NOPARTICLEINSERTION_H_
