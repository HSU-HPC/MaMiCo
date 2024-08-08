// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETKINETICENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETKINETICENERGYMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <cmath>
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class SetKineticEnergyMapping;
}
} // namespace coupling

/**
 *	@brief This class sets kinetic energy over several linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::SetKineticEnergyMapping {
public:
  /** obtains the old momentum over the region of interest. Besides,
   *  obtains the new momentum that shall be set.
   *	@param oldKineticEnergy
   *	@param newKineticEnergy
   *	@param numberParticles
   *	@param meanVelocity
   *	@param mdSolverInterface
   */
  SetKineticEnergyMapping(const double& oldKineticEnergy, const double& newKineticEnergy, const unsigned int& numberParticles,
                          const tarch::la::Vector<dim, double>& meanVelocity, coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _meanVelocity(meanVelocity),
        _correctionFactor(getCorrectionFactor(oldKineticEnergy, newKineticEnergy, numberParticles, meanVelocity)) {}

  /** Destructor */
  ~SetKineticEnergyMapping() {}

  /** empty function
   */
  void beginCellIteration() {}

  /** empty function
   */
  void endCellIteration() {}

  /** sets new velocity: still with same mean, but re-scale the deviation for
   *correct thermal energy
   *	@param cell
   */
  void handleCell(LinkedCell& cell) {
    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      coupling::interface::Molecule<dim>& wrapper(it->get());
      tarch::la::Vector<dim, double> velocity = wrapper.getVelocity();

      // set new velocity: still with same mean, but re-scale the deviation for
      // correct thermal energy
      velocity = _meanVelocity + _correctionFactor * (velocity - _meanVelocity);
      wrapper.setVelocity(velocity);

      it->next();
    }
    delete it;
  }

private:
  /** returns the correction factor between the old and new kinetic energy
   *	@param numberParticles
   *	@param momentum
   *	@return correctionFactor
   *	@remark no correction possible if the correction factor would tend to
   *infinity; I just hard-coded 1e-7 for this case
   */
  double getCorrectionFactor(const double& oldKineticEnergy, const double& newKineticEnergy, const unsigned int& numberParticles,
                             const tarch::la::Vector<dim, double>& meanVelocity) const {
    const double mass = _mdSolverInterface->getMoleculeMass();

    // no correction possible if the correction factor would tend to infinity; I
    // just hard-coded 1e-7 for this case
    if (oldKineticEnergy - 0.5 * mass * numberParticles * tarch::la::dot(meanVelocity, meanVelocity) < 1e-7) {
      return 1.0;
    }

    double correctionFactor = newKineticEnergy - 0.5 * mass * numberParticles * tarch::la::dot(meanVelocity, meanVelocity);
    correctionFactor = correctionFactor / (oldKineticEnergy - 0.5 * mass * numberParticles * tarch::la::dot(meanVelocity, meanVelocity));
    correctionFactor = sqrt(correctionFactor);
    return correctionFactor;
  }

  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  const tarch::la::Vector<dim, double> _meanVelocity;
  const double _correctionFactor;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETKINETICENERGYMAPPING_H_
