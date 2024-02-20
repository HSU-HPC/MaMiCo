// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_

#include "coupling/datastructures/CouplingCell.h"
#include "tarch/la/Vector.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class MomentumInsertion;
}

/** @brief used to manipulate the momentum/ velocity of the molecules contained
 * in a macroscopic cell.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 */
template <class LinkedCell, unsigned int dim> class coupling::MomentumInsertion {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface to the md solver*/
  MomentumInsertion(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface) : _mdSolverInterface(mdSolverInterface) {}
  /** @brief a simple destructor*/
  virtual ~MomentumInsertion() {}

  /** @brief returns the number of MD steps between subsequent momentum
   * insertions
   *  @returns the time step interval for momentum insertion */
  virtual unsigned int getTimeIntervalPerMomentumInsertion() const = 0;

  /** This method does not conserve the kinetic energy of the respective
   * macroscopic cell. To conserve the energy as well, see the description of
   * MomentumController on details how to do that.
   *  @brief inserts a fraction from the momentum of the macroscopic cell and
   * distributes is over all molecules.
   *  @param cell the macroscopic cell to insert the momentum
   *  @param fraction the fraction of momentum to use */
  virtual void insertMomentum(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell,
                              const unsigned int& currentMacroscopicCellIndex) const = 0;

protected:
  /** interface to the md solver */
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
};

#endif // _MOLECULARDYNAMICS_COUPLING_MOMENTUMINSERTION_H_
