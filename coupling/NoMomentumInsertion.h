// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_

#include "coupling/MomentumInsertion.h"
namespace coupling {
template <class LinkedCell, unsigned int dim> class NoMomentumInsertion;
}

/** @brief dummy class. Empty momentum insertion mechanism. Doesn't change
 * anything.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 */
template <class LinkedCell, unsigned int dim> class coupling::NoMomentumInsertion : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface to the md solver */
  NoMomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : MomentumInsertion<LinkedCell, dim>(mdSolverInterface) {}
  /** @brief a simple destructor */
  virtual ~NoMomentumInsertion() {}

  /** @brief returns the number of MD steps between subsequent momentum
   * insertions
   *  @todo We could set this to be zero? Since it does anything
   *  @returns the time step interval for momentum insertion */
  unsigned int getTimeIntervalPerMomentumInsertion() const override { return 1; }

  /** @brief a dummy function, which doesn't do anything
   * @param cell the coupling cell to not change */
  void insertMomentum(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 idx) const override {}

  void setInnerImposition(bool enable) override {}
};

#endif // _MOLECULARDYNAMICS_COUPLING_NOMOMENTUMINSERTION_H_
