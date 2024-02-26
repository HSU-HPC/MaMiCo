// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <iostream>
#include <list>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class ComputeKineticEnergyMapping;
}
} // namespace coupling

/**
 *	@brief This class computes the kinetic energy. inside a linked cell
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::ComputeKineticEnergyMapping {
public:
  /** Constructor
   *	@param mdSolverInterface
   */
  ComputeKineticEnergyMapping(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _kineticEnergy(0.0) {}

  /** Destructor */
  ~ComputeKineticEnergyMapping() {}

  /** sets the kinetic energy to zero, before the iteration process begins.
   */
  void beginCellIteration() { _kineticEnergy = 0.0; }

  /** computes the kinetic energy in a linked cell, by multiplying the specific
   * kinetic energy with the mass and 0.5.
   */
  void endCellIteration() { _kineticEnergy = 0.5 * _mdSolverInterface->getMoleculeMass() * _kineticEnergy; }

  /** sums up the specific kinetic energy of all molecules inside a linked cell.
   *	@param cell
   */
  void handleCell(LinkedCell& cell) {

    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      const coupling::interface::Molecule<dim>& wrapper(it->getConst());
      _kineticEnergy += tarch::la::dot(wrapper.getVelocity(), wrapper.getVelocity());

      it->next();
    }
    delete it;
  }

  /** returns the kinetic energy inside a linked cell
   *	@return _kineticEnergy
   */
  double getKineticEnergy() const { return _kineticEnergy; }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  double _kineticEnergy;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_
