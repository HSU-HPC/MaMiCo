// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_CELLMAPPINGS_PERTURBATEVELOCITYMAPPING_H_
#define _COUPLING_CELLMAPPINGS_PERTURBATEVELOCITYMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class PerturbateVelocityMapping;
}
} // namespace coupling

/** This class is used to superimpose a random perturbation to the mean velocity
 *in all directions. It is use for Dynamic MD cases to launch new MD instances
 *with independent molecular values (velocities)
 *	@brief This class is used to superimpose a random perturbation to the
 *mean velocity in all directions.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Niklas Wittmer
 */
template <class LinkedCell, unsigned int dim>
class coupling::cellmappings::PerturbateVelocityMapping {
public:
  /** Constructor
   *	@param mdSolverInterface
   *	@param velocity
   *	@param temperature
   */
  PerturbateVelocityMapping(
      coupling::interface::MDSolverInterface<LinkedCell, dim>
          *const mdSolverInterface,
      const tarch::la::Vector<dim, double> &velocity, const double &temperature)
      : _mdSolverInterface(mdSolverInterface),
        _molecularMass(mdSolverInterface->getMoleculeMass()),
        _kB(mdSolverInterface->getKB()), _temperature(temperature),
        _velocity(velocity), _sigma(_mdSolverInterface->getMoleculeSigma()) {}

  /** Destructor */
  ~PerturbateVelocityMapping() {}

  /** empty function
   */
  void beginCellIteration() {}
  /** empty function
   */
  void endCellIteration() {}

  /** superimposes a random perturbation to the mean velocity in all directions
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {

    double stdDeviation = std::sqrt(dim * _kB * _temperature / _molecularMass);

    tarch::la::Vector<dim, double> randomNumbers(0.0);

    coupling::interface::MoleculeIterator<LinkedCell, dim> *molecule =
        _mdSolverInterface->getMoleculeIterator(cell);
    molecule->begin();
    while (molecule->continueIteration()) {
      coupling::interface::Molecule<dim> &wrapper(molecule->get());
      randomNumbers[0] = tarch::utils::RandomNumberService::getInstance()
                             .getGaussianRandomNumber();
      for (unsigned int d = 1; d < dim; ++d) {
        randomNumbers[d] = tarch::utils::RandomNumberService::getInstance()
                               .getGaussianRandomNumber();
      }

      tarch::la::Vector<dim, double> mVelocity = wrapper.getVelocity();
      if (dim == 3) {
        mVelocity[0] =
            _velocity[0] +
            stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) *
                            std::cos(randomNumbers[2]));
        mVelocity[1] =
            _velocity[1] +
            stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) *
                            std::sin(randomNumbers[2]));
        mVelocity[2] =
            _velocity[2] +
            stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
      } else if (dim == 1) {
        mVelocity = _velocity + stdDeviation * randomNumbers;
      } else if (dim == 2) {
        mVelocity[0] =
            _velocity[0] +
            stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
        mVelocity[1] =
            _velocity[1] +
            stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]));
      }
      wrapper.setVelocity(mVelocity);

      molecule->next();
    }
    delete molecule;
  }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim>
      *const _mdSolverInterface;
  const double _molecularMass;
  const double _kB;
  const double _temperature;
  const tarch::la::Vector<dim, double> _velocity;
  const double _sigma;
};

#endif //_COUPLING_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_
