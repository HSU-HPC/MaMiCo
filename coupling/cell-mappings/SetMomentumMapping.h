// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETMOMENTUMMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETMOMENTUMMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class SetMomentumMapping;
}
} // namespace coupling

/**
 *	@brief This class sets a certain momentum over several linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::SetMomentumMapping {
public:
  /** obtains the old momentum over the region of interest. Besides,
   *  obtains the new momentum that shall be set and the number of particles
   *  contained in the macroscopic cell.
   *	@param oldMomentum
   *	@param newMomentum
   *	@param numberParticles
   *	@param mdSolverInterface
   */
  SetMomentumMapping(const tarch::la::Vector<dim, double>& oldMomentum, const tarch::la::Vector<dim, double>& newMomentum, const unsigned int& numberParticles,
                     coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _oldVelocity(getVelocity(numberParticles, oldMomentum)),
        _newVelocity(getVelocity(numberParticles, newMomentum)) {}

  /** Destructor */
  ~SetMomentumMapping() {}

  /** empty function
   */
  void beginCellIteration() {}

  /** empty function
   */
  void endCellIteration() {}

  /** applies a certain momentum over several linked cells, by steering the
   *velocity.
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      coupling::interface::Molecule<dim>& wrapper(it->get());
      tarch::la::Vector<dim, double> velocity = wrapper.getVelocity();
      velocity = velocity - _oldVelocity + _newVelocity;
      wrapper.setVelocity(velocity);

      it->next();
    }
    delete it;
  }

private:
  /** returns mean velocity og the cell
   *	@param numberParticles
   *	@param momentum
   *	@return mean velocity og the cell
   */
  tarch::la::Vector<dim, double> getVelocity(const unsigned int& numberParticles, const tarch::la::Vector<dim, double>& momentum) const {
    if (numberParticles != 0) {
      return (1.0 / (numberParticles * _mdSolverInterface->getMoleculeMass())) * momentum;
    } else {
      return tarch::la::Vector<dim, double>(0.0);
    }
  }

  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  tarch::la::Vector<dim, double> _oldVelocity;
  tarch::la::Vector<dim, double> _newVelocity;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETMOMENTUMMAPPING_H_
