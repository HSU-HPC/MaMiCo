// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREBINMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREBINMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class ComputeTemperatureBinMapping;
}
} // namespace coupling

/**
 *	@brief This class computes the momentum over certain linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Helene Wittenberg
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::ComputeTemperatureBinMapping {
public:
  /** Constructor
   *	@param mdSolverInterface
   */
  ComputeTemperatureBinMapping(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, tarch::la::Vector<dim, double>* velocity,
                               double* temperature)
      : _mdSolverInterface(mdSolverInterface), _velocity(velocity), _temperature(temperature) {}

  /** Destructor */
  ~ComputeTemperatureBinMapping() {}

  /** sets the mean velocity, momentum and the particle counter to zero, before the iteration process begins.
   */
  void beginCellIteration() {}

  /** computes the mean velocity, momentum in a linked cell,
   *	by dividing and multiplying the summation of the velocities computed in handleCell(...) over the number if particle and in particle mass respectively.
   */
  void endCellIteration() {}

  /** counts the molecules inside a linked cell and sums up the of the velocity of all particles inside the cell and saves it as momentum.
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      const coupling::interface::Molecule<dim>& wrapper(it->getConst());
      const int posZ = std::floor(wrapper.getPosition()[2] * 4);
      _temperature[posZ] += tarch::la::dot((wrapper.getVelocity() - _velocity[posZ]), (wrapper.getVelocity() - _velocity[posZ]));
      it->next();
    }
    delete it;
  }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim>* const _mdSolverInterface;
  tarch::la::Vector<dim, double>* _velocity;
  double* _temperature;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREBINMAPPING_H_
