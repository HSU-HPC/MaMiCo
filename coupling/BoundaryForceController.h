// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/datastructures/MacroscopicCell.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class BoundaryForceController;
}

/** controller for forces acting at open MD boundaries. There is an interface
 * method
 *  applyBoundaryForce which triggers potential boundary forcing in each
 * macroscopic
 *  cell that is located at the very outer MD boundary (first layer of non-ghost
 * macroscopic cells).
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::BoundaryForceController {
public:
  BoundaryForceController(coupling::interface::MDSolverInterface<
      LinkedCell, dim> *const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface) {}
  virtual ~BoundaryForceController() {}

  virtual void applyBoundaryForce(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,
                                                               dim> &cell,
      const unsigned int &currentLocalMacroscopicCellIndex) = 0;

  virtual double
  getPotentialEnergy(const tarch::la::Vector<dim, double> &position) const {
    return 0;
  }

  virtual tarch::la::Vector<dim, double>
  getForce(const tarch::la::Vector<dim, double> &position) const {
    return tarch::la::Vector<dim, double>(0.0);
  }

protected:
  coupling::interface::MDSolverInterface<LinkedCell,
                                         dim> *const _mdSolverInterface;
};
#endif // _MOLECULARDYNAMICS_COUPLING_BOUNDARYFORCECONTROLLER_H_
