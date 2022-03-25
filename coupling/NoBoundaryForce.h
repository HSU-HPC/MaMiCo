// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_
#define _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_

#include "coupling/BoundaryForceController.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class NoBoundaryForce;
}

/** @brief dummy implementation, applying no boundary force.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 */
template <class LinkedCell, unsigned int dim> class coupling::NoBoundaryForce : public coupling::BoundaryForceController<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface to the md solver*/
  NoBoundaryForce(coupling::interface::MDSolverInterface<LinkedCell, dim> *mdSolverInterface)
      : coupling::BoundaryForceController<LinkedCell, dim>(mdSolverInterface) {}
  /** @brief a simple destructor*/
  virtual ~NoBoundaryForce() {}

  /** @brief dummy function, doesn't do anything
   *  @param cell the macroscopic cell to apply no force
   *  @param currentLocalMacroscopicCellIndex the linearised local index of the
   * macroscopic cell*/
  virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                                  const unsigned int &currentLocalMacroscopicCellIndex) {}
};
#endif // _MOLECULARDYNAMICS_COUPLING_NOBOUNDARYFORCE_H_
