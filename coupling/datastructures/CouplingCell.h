// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATA_STRUCTURES_COUPLINGCELL_H_
#define _MOLECULARDYNAMICS_COUPLING_DATA_STRUCTURES_COUPLINGCELL_H_

#include "tarch/la/Vector.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace coupling {
namespace datastructures {

template <unsigned int dim> class CouplingCell;
} // namespace datastructures
} // namespace coupling

/** describes a quadratic/ cubic coupling cell filled with fluid (no linked
 *cells). Base class for the class
 *coupling::datastructures::CouplingCellWithLinkedCells
 *	@brief defines the cell type with cell-averaged quantities only (no
 *linked cells).
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::datastructures::CouplingCell {
public:
  /** Constructor: initialises the coupling cell with zero values.
   */
  CouplingCell()
      : _microscopicMass(0.0), _microscopicMomentum(0.0), _macroscopicMass(0.0), _macroscopicMomentum(0.0), _potentialEnergy(0.0), _temperature(0.0),
        _currentVelocity(0.0) {}

  /** Destructor */
  virtual ~CouplingCell() {}

  /** sets the microscopic mass
   * @param mass Mass*/
  void setMicroscopicMass(const double& mass) { _microscopicMass = mass; }
  /** returns the microscopic mass
   * @returns _microscopicMass Mass*/
  const double& getMicroscopicMass() const { return _microscopicMass; }
  /** sets the microscopic moments
   * @param momentum Momentum*/
  void setMicroscopicMomentum(const tarch::la::Vector<dim, double>& momentum) { _microscopicMomentum = momentum; }
  /** returns the microscopic moments
   * @returns _microscopicMomentum Momentum*/
  const tarch::la::Vector<dim, double>& getMicroscopicMomentum() const { return _microscopicMomentum; }

  /** sets the microscopic mass
   * @param mass Mass*/
  void setMacroscopicMass(const double& mass) { _macroscopicMass = mass; }
  /** returns the microscopic mass
   * @returns _microscopicMass Mass*/
  const double& getMacroscopicMass() const { return _macroscopicMass; }
  /** sets the microscopic moments
   * @param momentum Momentum*/
  void setMacroscopicMomentum(const tarch::la::Vector<dim, double>& momentum) { _macroscopicMomentum = momentum; }
  /** returns the microscopic moments
   * @returns _microscopicMomentum Momentum*/
  const tarch::la::Vector<dim, double>& getMacroscopicMomentum() const { return _macroscopicMomentum; }

  /** returns the mean potential energy over the coupling cell
   * @returns _potentialEnergy potential energy*/
  const double& getPotentialEnergy() const { return _potentialEnergy; }
  /** sets the mean potential energy over the coupling cell
   * @param potentialEnergy potential energy*/
  void setPotentialEnergy(const double& potentialEnergy) { _potentialEnergy = potentialEnergy; }

  /** adds a certain amount to the microscopic mass
   * @param mass Mass to be added to the microscopic mass*/
  void addMicroscopicMass(const double& mass) { _microscopicMass += mass; }
  /** adds a certain amount to the microscopic moments
   * @param momentum Momentum to be added to the microscopic moments*/
  void addMicroscopicMomentum(const tarch::la::Vector<dim, double>& momentum) { _microscopicMomentum += momentum; }

  /** adds a certain amount to the macroscopic mass
   * @param mass Mass to be added to the macroscopic mass*/
  void addMacroscopicMass(const double& mass) { _macroscopicMass += mass; }
  /** adds a certain amount to the macroscopic moments
   * @param momentum Momentum to be added to the macroscopic moments*/
  void addMacroscopicMomentum(const tarch::la::Vector<dim, double>& momentum) { _macroscopicMomentum += momentum; }

  /** sets current velocity (sampled right before any distributeX(...) call)
   * @param velocity velocity*/
  void setCurrentVelocity(const tarch::la::Vector<dim, double>& velocity) { _currentVelocity = velocity; }
  /** returns current velocity (sampled right before any distributeX(...) call)
   * @return _currentVelocity velocity*/
  const tarch::la::Vector<dim, double>& getCurrentVelocity() const { return _currentVelocity; }

  /** sets the temperature
   * @param temperature temperature*/
  void setTemperature(const double& temperature) { _temperature = temperature; }
  /** returns the temperature
   * @return _temperature temperature*/
  const double& getTemperature() const { return _temperature; }

  /** buffers for macroscopic mass that need to be transferred
   *  from the macroscopic to the microscopic simulation */
  double _microscopicMass;
  /** buffers for macroscopic quantities of momentum that need to be transferred
   *  from the macroscopic to the microscopic simulation */
  tarch::la::Vector<dim, double> _microscopicMomentum;

  /** buffers for macroscopic mass that need to be returned
   *  from the microscopic to the macroscopic simulation */
  double _macroscopicMass;
  /** buffers for macroscopic quantities of momentum that need to be returned
   *  from the microscopic to the macroscopic simulation */
  tarch::la::Vector<dim, double> _macroscopicMomentum;

  /** holds the mean potential energy within the coupling cell. This value is
   * needed as a reference potential energy value for the USHER scheme.
   */
  double _potentialEnergy;

  /** temperature within this cell. Needed for the thermostat. */
  double _temperature;

  /** buffer for current mean velocity in the cell. */
  tarch::la::Vector<dim, double> _currentVelocity;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATA_STRUCTURES_COUPLINGCELL_H_
