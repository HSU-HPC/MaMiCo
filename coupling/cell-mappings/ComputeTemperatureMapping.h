// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class ComputeTemperatureMapping;
}
} // namespace coupling

/**
 *	@brief This class computes the temperature in a certain (macroscopic) cell.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::ComputeTemperatureMapping {
public:
  /** Constructor
   *	@param meanVelocity
   *	@param mdSolverInterface
   */
  ComputeTemperatureMapping(const tarch::la::Vector<dim, double> &meanVelocity,
                            coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _meanVelocity(meanVelocity), _temperature(0.0), _particleCounter(0) {}

  /** Destructor */
  ~ComputeTemperatureMapping() {}

  /** sets the temperature and the particle counter to zero, before the
iteration process begins.
         */
  void beginCellIteration() {
    _temperature = 0.0;
    _particleCounter = 0;
  }

  /** computes the temperature in a linked cell, based on equilibrium
statistical mechanics.
         */
  void endCellIteration() {
    _temperature = _temperature * _mdSolverInterface->getMoleculeMass();
    if (_particleCounter != 0) {
      _temperature = _temperature / (dim * _mdSolverInterface->getKB() * _particleCounter);
    }
  }

  /** sums up the velocity fluctuation (from the mean flow velocity) of all
 particles
         *	@param cell
         *	@param cellIndex
         */
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      const coupling::interface::Molecule<dim> &wrapper(it->getConst());
      _temperature += tarch::la::dot((wrapper.getVelocity() - _meanVelocity), (wrapper.getVelocity() - _meanVelocity));
      _particleCounter++;

      it->next();
    }
    delete it;
  }

  /** returns the temperature inside a linked cell
   *	@return _temperature
   */
  double getTemperature() const { return _temperature; }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim> *const _mdSolverInterface;
  const tarch::la::Vector<dim, double> _meanVelocity;
  double _temperature;
  unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
