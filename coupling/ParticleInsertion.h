// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_

#include "coupling/BoundaryForceController.h"
#include <cstdlib>
#include <iostream>

namespace coupling {
template <class LinkedCell, unsigned int dim> class ParticleInsertion;
}

/** @brief interface for particle insertion/removal on macroscopic cell basis.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 */
template <class LinkedCell, unsigned int dim>
class coupling::ParticleInsertion {
public:
  /** this state is returned by the insertDeleteMass() function. To tell the
   * program, which action happened through the application of the function*/
  enum Action {
    NoAction = 0,  ///< No action was taken / Nothing changed
    Insertion = 1, ///< A molecule was inserted
    Deletion = 2   ///< A molecule was deleted
  };

  /**
   *  @brief adds or removes particles to the macroscopic cell, simulates a mass
   * flow.
   *  @param cell the macroscopic cell to insert or delete mass
   *  @param macroscopicCellPosition position of the cell given as a vector
   *  @param macroscopicCellSize size of the macroscopic cell given as a vector
   *  @param meanVelocity the mean velocity in the cell as a vector
   *  @param temperature the current temperature in the cell
   *  @param boundaryForceController an instance of the boundaryForceController
   * for the simulation
   *  @returns the type of coupling::ParticleInsertion::Action that was applied
   */
  virtual typename coupling::ParticleInsertion<LinkedCell, dim>::Action
  insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>
          &cell,
      const tarch::la::Vector<dim, double> &macroscopicCellPosition,
      const tarch::la::Vector<dim, double> &macroscopicCellSize,
      const tarch::la::Vector<dim, double> &meanVelocity,
      const double &temperature,
      const coupling::BoundaryForceController<LinkedCell, dim>
          &boundaryForceController) = 0;

  /** returns true, if the particle insertion requires information on the
   * potential energy landscape. The USHER algorithm requires it. Other
   * algorithms may not require it; the trivial
   * NoParticleInsertion-implementation which does not do anything obviously
   * does not require a valid potential energy landscape.
   *  @brief returns true, if a potential energy landscape is needed for the
   * insertion/removal.
   *  @returns a bool that indicates if a potential energy landscape is
   * neccessary (true) or not (false) */
  virtual bool requiresPotentialEnergyLandscape() = 0;

  /** @brief a simple constructor
   *  @param insertDeleteMassAtTimestep interval of time steps for the
   * insertion/removal of particles*/
  ParticleInsertion(unsigned int insertDeleteMassEveryTimestep)
      : _insertDeleteMassEveryTimestep(insertDeleteMassEveryTimestep) {
    if (_insertDeleteMassEveryTimestep == 0) {
      std::cout << "ERROR ParticleInsertion::ParticleInsertion(..): "
                   "_insertDeleteMassEveryTimestep=0!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** @brief a simple destructor*/
  virtual ~ParticleInsertion() {}

  /** @brief returns true if mass needs to be inserted or removed in a time step
   * t
   *  @param t time step to check for the insertion/removal of mass
   *  @returns a bool that indicates if mass insertion/deletion shall be applied
   * in this time step (true) or not (false) */
  bool insertDeleteMassAtTimestep(unsigned int t) const {
    return (t % _insertDeleteMassEveryTimestep == 0);
  }

private:
  /** interval of time steps for the insertion/removal of particles */
  const unsigned int _insertDeleteMassEveryTimestep;
};
#endif // _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
