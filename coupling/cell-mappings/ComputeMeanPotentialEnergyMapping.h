// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_

#include <iostream>
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/BoundaryForceController.h"

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim>
class ComputeMeanPotentialEnergyMapping;
}
}

/**
 *	@brief This class computes the mean potential energy over this macroscopic
 * cell.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::cellmappings::ComputeMeanPotentialEnergyMapping {
public:
  /** Constructor
	 *	@param mdSolverInterface
	 *	@param boundaryForceController
	 */
  ComputeMeanPotentialEnergyMapping(
      coupling::interface::MDSolverInterface<LinkedCell,
                                             dim> *const mdSolverInterface,
      const coupling::BoundaryForceController<LinkedCell, dim> &
          boundaryForceController)
      : _mdSolverInterface(mdSolverInterface), _meanPotentialEnergy(0.0),
        _particleCounter(0), _boundaryForceController(boundaryForceController) {
  }

  /** Destructor */
  ~ComputeMeanPotentialEnergyMapping() {}

  /** sets the mean potential energy and the particle counter to zero, before
the iteration process begins.
	 */
  void beginCellIteration() {
    _meanPotentialEnergy = 0.0;
    _particleCounter = 0;
  }

  /** computes the mean potential energy in a linked cell, by dividing the
summation of the mean potential energy of all particles inside the cell over the
number of particles.
	 */
  void endCellIteration() {
    if (_particleCounter != 0) {
      _meanPotentialEnergy = _meanPotentialEnergy / _particleCounter;
    }
  }

  /** counts the molecules inside a linked cell and sums up the of the mean
 potential energy of all particles inside the cell.
 	 *	@param cell
 	 *	@param cellIndex
 	 */
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim> *it =
        _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      const coupling::interface::Molecule<dim> &wrapper(it->getConst());
      _meanPotentialEnergy += wrapper.getPotentialEnergy();
      _meanPotentialEnergy +=
          _boundaryForceController.getPotentialEnergy(wrapper.getPosition());
      _particleCounter++;

      it->next();
    }
    delete it;
  }

  /** returns the mean potential energy inside a linked cell
	 *	@return _meanPotentialEnergy
	 */
  double getPotentialEnergy() const { return _meanPotentialEnergy; }

private:
  coupling::interface::MDSolverInterface<LinkedCell,
                                         dim> *const _mdSolverInterface;
  double _meanPotentialEnergy;
  unsigned int _particleCounter;
  const coupling::BoundaryForceController<LinkedCell, dim> &
      _boundaryForceController;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_
