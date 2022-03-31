// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_

#include "tarch/la/Vector.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/SetMomentumMapping.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class MomentumController;
}

/** controles the momentum in a macroscopic cell. This class can compute the
 * momentum and set
 *  a certain value for the momentum.
 *
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::MomentumController {
public:
  MomentumController(coupling::interface::MDSolverInterface<
      LinkedCell, dim> *const mdSolverInterface)
      : _computeMassMapping(mdSolverInterface),
        _computeMomentumMapping(mdSolverInterface),
        _mdSolverInterface(mdSolverInterface) {}
  ~MomentumController() {}

  /** computes and returns the momentum and mean velocity of a macroscopic cell
   */
  void computeMomentumAndMeanVelocity(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,
                                                               dim> &cell,
      tarch::la::Vector<dim, double> &momentum,
      tarch::la::Vector<dim, double> &meanVelocity) {
    cell.iterateConstCells(_computeMomentumMapping);
    momentum = _computeMomentumMapping.getMomentum();
    meanVelocity = _computeMomentumMapping.getMeanVelocity();
  }

  /** computes and returns the momentum of a macroscopic cell */
  void computeMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<
                           LinkedCell, dim> &cell,
                       tarch::la::Vector<dim, double> &momentum) {
    cell.iterateConstCells(_computeMomentumMapping);
    momentum = _computeMomentumMapping.getMomentum();
  }

  /** computes and returns the mean velocity of a macroscopic cell */
  void
  computeMeanVelocity(coupling::datastructures::MacroscopicCellWithLinkedCells<
                          LinkedCell, dim> &cell,
                      tarch::la::Vector<dim, double> &meanVelocity) {
    cell.iterateConstCells(_computeMomentumMapping);
    meanVelocity = _computeMomentumMapping.getMeanVelocity();
  }

  /** sets the momentum of a macroscopic cell to newMomentum. This cell does not
   * consider the kinetic energy
   *  of the system, so the kinetic energy will be different after executing
   * this method!
   *  In order to also retain the kinetic energy, one might compute the kinetic
   * energy before calling this method
   *  (using the KineticEnergyController) and set the old kinetic energy
   * afterwards (again using the KineticEnergyController).
   *  This will then set the right kinetic energy again, conserving both mass
   * and momentum over this cell.
   */
  void setMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<
                       LinkedCell, dim> &cell,
                   const tarch::la::Vector<dim, double> &newMomentum) {
    tarch::la::Vector<dim, double> currentMomentum(0.0);
    computeMomentum(cell, currentMomentum);
    cell.iterateConstCells(_computeMassMapping);
    unsigned int numberParticles = _computeMassMapping.getNumberOfParticles();
    coupling::cellmappings::SetMomentumMapping<LinkedCell, dim>
        setMomentumMapping(currentMomentum, newMomentum, numberParticles,
                           _mdSolverInterface);
    cell.iterateCells(setMomentumMapping);
  }

private:
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim>
      _computeMassMapping;
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
      _computeMomentumMapping;
  coupling::interface::MDSolverInterface<LinkedCell,
                                         dim> *const _mdSolverInterface;
};
#endif // _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_
