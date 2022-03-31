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

/** applies the interpolating boundary force from Zhou et al. in every boundary
 * cell. For details on the forcing,
 *  check out the descriptions in cell-mappings/ZhouBoundaryForce.
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::ZhouBoundaryForceController : public coupling::BoundaryForceController<LinkedCell, dim> {
public:
  ZhouBoundaryForceController(const double &density, const double &temperature, const tarch::la::Vector<2 * dim, bool> &boundary,
                              coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface)
      : coupling::BoundaryForceController<LinkedCell, dim>(mdSolverInterface), _density(density), _temperature(temperature), _boundary(boundary),
        _zhouBoundaryForce(density, temperature, mdSolverInterface->getMoleculeEpsilon(), mdSolverInterface->getMoleculeSigma(), boundary,
                           mdSolverInterface->getGlobalMDDomainOffset(), mdSolverInterface->getGlobalMDDomainSize(), mdSolverInterface) {}

  virtual ~ZhouBoundaryForceController() {}

  virtual void applyBoundaryForce(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                                  const unsigned int &currentLocalMacroscopicCellIndex) {
    cell.iterateCells(_zhouBoundaryForce);
  }

  virtual double getPotentialEnergy(const tarch::la::Vector<dim, double> &position) const { return _zhouBoundaryForce.getPotentialEnergy(position); }

  virtual tarch::la::Vector<dim, double> getForce(const tarch::la::Vector<dim, double> &position) const {
    return _zhouBoundaryForce.getBoundaryForces(position);
  }

private:
  const double _density;
  const double _temperature;
  const tarch::la::Vector<2 * dim, bool> _boundary;
  coupling::cellmappings::ZhouBoundaryForce<LinkedCell, dim> _zhouBoundaryForce;
};
#endif // _MOLECULARDYNAMICS_COUPLING_ZHOUBOUNDARYFORCECONTROLLER_H_
