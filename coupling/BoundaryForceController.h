// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_

#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class BoundaryForceController;
}

/** There is an interface method applyBoundaryForce which triggers potential
 * boundary forcing in each macroscopic cell that is located at the very outer
 * MD boundary (first layer of non-ghost macroscopic cells).
 *  @brief controller for forces acting at open MD boundaries
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim the integer dim refers to the spacial dimension of the
 * simulation, can be 1, 2, or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::BoundaryForceController {
public:
  /**@brief A simple constructor
   *  @param mdSolverInterface interface to the molecular dynamics solver*/
  BoundaryForceController(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface) : _mdSolverInterface(mdSolverInterface) {}
  /**@brief A simple destructor*/
  virtual ~BoundaryForceController() {}

  /** iterates over all linked cells of the given macroscopic cell and applies
   * the cellmapping for the boundary force
   * @brief applies the boundary force on a boundary cell
   * @param cell the macroscopic boundary cell to apply the boundary force
   * @param currentLocalMacroscopicCellIndex the index of the macroscopic cell
   */
  virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell,
                                  const unsigned int& currentLocalMacroscopicCellIndex) = 0;

  /** @brief calculates the potential energy for a given position
   *  @param position the position for which the potential energy will be
   * calculated
   *  @returns the potential energy for the given position  */
  virtual double getPotentialEnergy(const tarch::la::Vector<dim, double>& position) const { return 0; }

  /** @brief calculates the boundary force for the given particle position
   *  @param position particle position for the force calculation
   *  @returns the force for the given position */
  virtual tarch::la::Vector<dim, double> getForce(const tarch::la::Vector<dim, double>& position) const { return tarch::la::Vector<dim, double>(0.0); }

protected:
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface; ///< interface of the molecular dynamics solver
};
#endif // _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_
