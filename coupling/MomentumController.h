// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/SetMomentumMapping.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "tarch/la/Vector.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class MomentumController;
}

/** This class can compute the momentum and set a certain value for the
 * momentum.
 *  @brief controls the momentum in a macroscopic cell.
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3 */
template <class LinkedCell, unsigned int dim> class coupling::MomentumController {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface to the md solver */
  MomentumController(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : _computeMassMapping(mdSolverInterface), _computeMomentumMapping(mdSolverInterface), _mdSolverInterface(mdSolverInterface) {}

  /** a simple destructor */
  ~MomentumController() {}

  /** @brief computes and returns the momentum and mean velocity of a
   * macroscopic cell
   *  @param cell the macroscopic cell, for which the values shall be calculated
   *  @param momentum vector to which the momentum of the cell will be written
   * to
   *  @param meanVelocity vector to which the mean velocity of the cell will be
   * written to*/
  void computeMomentumAndMeanVelocity(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell, tarch::la::Vector<dim, double>& momentum,
                                      tarch::la::Vector<dim, double>& meanVelocity) {
    cell.iterateConstCells(_computeMomentumMapping);
    momentum = _computeMomentumMapping.getMomentum();
    meanVelocity = _computeMomentumMapping.getMeanVelocity();
  }

  /** @brief computes and returns the momentum of a macroscopic cell
   * @param cell the macroscopic cell, for which the values shall be calculated
   *  @param momentum vector to which the momentum of the cell will be written
   * to */
  void computeMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell, tarch::la::Vector<dim, double>& momentum) {
    cell.iterateConstCells(_computeMomentumMapping);
    momentum = _computeMomentumMapping.getMomentum();
  }

  /** computes and returns the mean velocity of a macroscopic cell
   *  @param cell the macroscopic cell, for which the values shall be calculated
   *  @param meanVelocity vector to which the mean velocity of the cell will be
   * written to*/
  void computeMeanVelocity(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell, tarch::la::Vector<dim, double>& meanVelocity) {
    cell.iterateConstCells(_computeMomentumMapping);
    meanVelocity = _computeMomentumMapping.getMeanVelocity();
  }

  /** This cell does not consider the kinetic energy of the system, so the
   * kinetic energy will be different after executing this method! In order to
   * also retain the kinetic energy, one might compute the kinetic energy before
   * calling this method (using the KineticEnergyController) and set the old
   * kinetic energy afterwards (again using the KineticEnergyController). This
   * will then set the right kinetic energy again, conserving both mass and
   * momentum over this cell.
   *  @brief sets the momentum of a macroscopic cell to the input value.
   *  @param cell the macroscopic cell in which the momentum shall be changed
   *  @param newMomentum the value to which the momentum will be set */
  void setMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim>& cell, const tarch::la::Vector<dim, double>& newMomentum) {
    tarch::la::Vector<dim, double> currentMomentum(0.0);
    computeMomentum(cell, currentMomentum);
    cell.iterateConstCells(_computeMassMapping);
    unsigned int numberParticles = _computeMassMapping.getNumberOfParticles();
    coupling::cellmappings::SetMomentumMapping<LinkedCell, dim> setMomentumMapping(currentMomentum, newMomentum, numberParticles, _mdSolverInterface);
    cell.iterateCells(setMomentumMapping);
  }

private:
  /** instance of the ComputeMassMapping */
  coupling::cellmappings::ComputeMassMapping<LinkedCell, dim> _computeMassMapping;
  /** instance of the ComputeMomentumMapping*/
  coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> _computeMomentumMapping;
  /** interface to the md solver*/
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
};
#endif // _MOLECULARDYNAMICS_COUPLING_MOMENTUMCONTROLLER_H_
