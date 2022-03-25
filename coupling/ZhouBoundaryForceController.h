// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_

#include "coupling/BoundaryForceController.h"
#include "coupling/cell-mappings/ZhouBoundaryForce.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class ZhouBoundaryForceController;
}

/** For details on the forcing, check out the descriptions in
 * cell-mappings/ZhouBoundaryForce.
 *  @brief applies the boundary force from Zhou et al. in boundary cell.
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim the integer dim refers to the spacial dimension of the
 * simulation, can be 1, 2, or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::ZhouBoundaryForceController : public coupling::BoundaryForceController<LinkedCell, dim> {
public:
  /** @brief A simple constructor
   * @param density the mean density for the setup
   * @param temperature the mean temperature for the setup
   * @param boundary a vector of booleans indicating at which boundaries the
   * force shall be applied
   * @param mdSolverInterface a reference to the interface of the MD solver*/
  ZhouBoundaryForceController(const double &density, const double &temperature, const tarch::la::Vector<2 * dim, bool> &boundary,
                              coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface)
      : coupling::BoundaryForceController<LinkedCell, dim>(mdSolverInterface), _density(density), _temperature(temperature), _boundary(boundary),
        _zhouBoundaryForce(density, temperature, mdSolverInterface->getMoleculeEpsilon(), mdSolverInterface->getMoleculeSigma(), boundary,
                           mdSolverInterface->getGlobalMDDomainOffset(), mdSolverInterface->getGlobalMDDomainSize(), mdSolverInterface) {}

  /** @brief Destructor*/
  virtual ~ZhouBoundaryForceController() {}

  /** iterates over all linked cells of the given macroscopic cell and applies
   * the cellmapping for the Zhou boundary force
   * @brief applies the Zhou boundary force on a boundary cell
   * @param cell the macroscopic boundary cell to apply the boundary force
   * @param currentLocalMacroscopicCellIndex the index of the macroscopic cell
   */
  virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                                  const unsigned int &currentLocalMacroscopicCellIndex) {
    cell.iterateCells(_zhouBoundaryForce);
  }

  /** @brief calculates the potential energy for a given position
   *  @param position the position for which the potential energy will be
   * calculated
   *  @returns the potential energy at the given position */
  virtual double getPotentialEnergy(const tarch::la::Vector<dim, double> &position) const { return _zhouBoundaryForce.getPotentialEnergy(position); }

  /** @brief calculates the boundary force for the given particle position
   *  @param position particle position for the force calculation
   *  @returns the boundary force at the given position */
  virtual tarch::la::Vector<dim, double> getForce(const tarch::la::Vector<dim, double> &position) const {
    return _zhouBoundaryForce.getBoundaryForces(position);
  }

private:
  const double _density;     ///< the mean density
  const double _temperature; ///< the mean temperature
  /** A boolean vector which indicates where to apply the boundary force.
   * Per dimension it has two enteries. They refer to the boundaries in the
   * following order: left (small x), right (high x), front (small y), back
   * (high y), bottom (small z), top (high z).
   * @brief indicates at which boundaries to apply the force  */
  const tarch::la::Vector<2 * dim, bool> _boundary;
  coupling::cellmappings::ZhouBoundaryForce<LinkedCell, dim> _zhouBoundaryForce; ///< the cell mapping for the application of the force
};
#endif // _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_
