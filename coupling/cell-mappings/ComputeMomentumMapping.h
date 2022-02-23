// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class ComputeMomentumMapping;
}
} // namespace coupling

/**
 *	@brief This class computes the momentum over certain linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::cellmappings::ComputeMomentumMapping {
public:
  /** Constructor
   *	@param mdSolverInterface
   */
  ComputeMomentumMapping(coupling::interface::MDSolverInterface<LinkedCell, dim>
                             *const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _momentum(0.0),
        _meanVelocity(0.0), _particleCounter(0) {}

  /** Destructor */
  ~ComputeMomentumMapping() {}

  /** sets the mean velocity, momentum and the particle counter to zero, before
   * the iteration process begins.
   */
  void beginCellIteration() {
    _momentum = tarch::la::Vector<dim, double>(0.0);
    _meanVelocity = tarch::la::Vector<dim, double>(0.0);
    _particleCounter = 0;
  }

  /** computes the mean velocity, momentum in a linked cell,
   *	by dividing and multiplying the summation of the velocities computed in
   *handleCell(...) over the number if particle and in particle mass
   *respectively.
   */
  void endCellIteration() {
    if (_particleCounter != 0) {
      _meanVelocity = (1.0 / ((double)_particleCounter)) * _momentum;
      _momentum = _mdSolverInterface->getMoleculeMass() * _momentum;
    }
  }

  /** counts the molecules inside a linked cell and sums up the of the velocity
   *of all particles inside the cell and saves it as momentum.
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim> *it =
        _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      const coupling::interface::Molecule<dim> &wrapper(it->getConst());
      _momentum += wrapper.getVelocity();
      _particleCounter++;

      it->next();
    }
    delete it;
  }

  /** returns the momentum inside a linked cell
   *	@return _momentum
   */
  tarch::la::Vector<dim, double> getMomentum() const { return _momentum; }

  /** returns the mean velocity inside a linked cell
   *	@return _meanVelocity
   */
  tarch::la::Vector<dim, double> getMeanVelocity() const {
    return _meanVelocity;
  }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim>
      *const _mdSolverInterface;
  tarch::la::Vector<dim, double> _momentum;
  tarch::la::Vector<dim, double> _meanVelocity;
  unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_
