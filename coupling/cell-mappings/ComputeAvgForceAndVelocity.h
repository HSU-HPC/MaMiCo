// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_

#include "coupling/interface/MDSolverInterface.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class ComputeAvgForceAndVelocity;
}
} // namespace coupling

/** This class sums up all force and velocity vectors and counts molecules
 *inside a linked cell. Afterwards, the average force/velocity contribution is
 *computed.
 *	@brief This class sums up all force/velocity vectors and counts
 *molecules inside a linked cell
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::ComputeAvgForceAndVelocity {
public:
  /** Constructor
   *	@param mdSolverInterface
   */
  ComputeAvgForceAndVelocity(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _force(0.0), _velocity(0.0), _particleCounter(0) {}

  /** Destructor */
  ~ComputeAvgForceAndVelocity() {}

  /** sets force, velocity and moluce counter to zero, before the iteration
   * process begins.
   */
  void beginCellIteration() {
    _particleCounter = 0;
    _force = tarch::la::Vector<dim, double>(0.0);
    _velocity = tarch::la::Vector<dim, double>(0.0);
  }

  /** the average force and velocity contribution incide a linked cell are
   *computed, by dividing the summation calculated in endCellIteration() over
   *number of the particles inside the cell.
   */
  void endCellIteration() {
    if (_particleCounter != 0) {
      _force = (1.0 / _particleCounter) * _force;
      _velocity = (1.0 / _particleCounter) * _velocity;
    }
  }

  /** sums up all force and velocity vectors and counts molecules inside a
   *linked cell.
   *	@param cell
   */
  void handleCell(LinkedCell& cell) {
    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      coupling::interface::Molecule<dim>& wrapper(it->get());
      _particleCounter++;
      _force = _force + wrapper.getForce();
      _velocity = _velocity + wrapper.getVelocity();
      it->next();
    }
    delete it;
  }

  /** returns the force vectors inside a linked cell
   *	@return _force
   */
  tarch::la::Vector<dim, double> getAvgForce() const { return _force; }

  /** returns the velocity vectors inside a linked cell
   *	@return _velocity
   */
  tarch::la::Vector<dim, double> getAvgVelocity() const { return _velocity; }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  tarch::la::Vector<dim, double> _force;
  tarch::la::Vector<dim, double> _velocity;
  unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_
