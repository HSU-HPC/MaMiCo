// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETTEMPERATUREMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETTEMPERATUREMAPPING_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include <cmath>
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class SetTemperatureMapping;
}
} // namespace coupling

/**
 *	@brief This class sets a certain temperature over several linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::SetTemperatureMapping {
public:
  /** Constructor: obtains the old temperature over the region of interest.
Besides,
     *  obtains the new temperature.
         *	@param oldTemperature
         *	@param newTemperature
         *	@param meanVelocity
         *	@param mdSolverInterface
     */
  SetTemperatureMapping(const double &oldTemperature, const double &newTemperature, const tarch::la::Vector<dim, double> &meanVelocity,
                        coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface), _factor(getScalingFactor(oldTemperature, newTemperature)), _meanVelocity(meanVelocity) {}

  /** Destructor */
  ~SetTemperatureMapping() {}

  /** empty function
   */
  void beginCellIteration() {}

  /** empty function
   */
  void endCellIteration() {}

  /** applies a certain temperature over several linked cells, by changing the
velocity (velocity fluctuation).
         *	@param cell
         *	@param cellIndex
         */
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      coupling::interface::Molecule<dim> &wrapper(it->get());
      tarch::la::Vector<dim, double> velocity = wrapper.getVelocity();
      // scale velocity
      velocity = _meanVelocity + _factor * (velocity - _meanVelocity);
      wrapper.setVelocity(velocity);
      it->next();
    }
    delete it;
  }

private:
  /** calculates the scaling factor between the old and new temperatures.
         *	@param oldTemperature
         *	@param newTemperature
         *	@return newTemperature
         *	@remark only allow re-scaling if the original temperature was
not zero (can happen in empty cells that just became populated)
         */
  double getScalingFactor(const double &oldTemperature, const double &newTemperature) const {

    if (oldTemperature != 0.0) {
      return sqrt(newTemperature / oldTemperature);
    } else {
      return 1.0;
    }
  }

  coupling::interface::MDSolverInterface<LinkedCell, dim> *const _mdSolverInterface;
  const double _factor;
  const tarch::la::Vector<dim, double> _meanVelocity;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_SETTEMPERATUREMAPPING_H_
