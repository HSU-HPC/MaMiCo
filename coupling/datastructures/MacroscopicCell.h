// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELL_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELL_H_

#include "tarch/la/Vector.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace coupling {
namespace datastructures {

template <unsigned int dim> class MacroscopicCell;
template <class LinkedCell, unsigned int dim> class MacroscopicCellWithLinkedCells;
} // namespace datastructures
} // namespace coupling

/** describes a quadratic/ cubic macroscopic cell filled with fluid (no linked
 * cells). Base class for the class
 * coupling::datastructures::MacroscopicCellWithLinkedCells
 *	@brief defines the cell type with cell-averaged quantities only (no
 *linked cells).
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::datastructures::MacroscopicCell {
public:
  /** Constructor: initialises the macroscopic cell with zero values.
   */
  MacroscopicCell()
      : _microscopicMass(0.0), _microscopicMomentum(0.0), _macroscopicMass(0.0), _macroscopicMomentum(0.0), _potentialEnergy(0.0), _temperature(0.0),
        _currentVelocity(0.0) {}

  /** Destructor */
  virtual ~MacroscopicCell() {}

  /** sets the microscopic mass
   * @param mass Mass*/
  void setMicroscopicMass(const double &mass) { _microscopicMass = mass; }
  /** returns the microscopic mass
   * @returns _microscopicMass Mass*/
  const double &getMicroscopicMass() const { return _microscopicMass; }
  /** sets the microscopic moments
   * @param momentum Momentum*/
  void setMicroscopicMomentum(const tarch::la::Vector<dim, double> &momentum) { _microscopicMomentum = momentum; }
  /** returns the microscopic moments
   * @returns _microscopicMomentum Momentum*/
  const tarch::la::Vector<dim, double> &getMicroscopicMomentum() const { return _microscopicMomentum; }

  /** sets the microscopic mass
   * @param mass Mass*/
  void setMacroscopicMass(const double &mass) { _macroscopicMass = mass; }
  /** returns the microscopic mass
   * @returns _microscopicMass Mass*/
  const double &getMacroscopicMass() const { return _macroscopicMass; }
  /** sets the microscopic moments
   * @param momentum Momentum*/
  void setMacroscopicMomentum(const tarch::la::Vector<dim, double> &momentum) { _macroscopicMomentum = momentum; }
  /** returns the microscopic moments
   * @returns _microscopicMomentum Momentum*/
  const tarch::la::Vector<dim, double> &getMacroscopicMomentum() const { return _macroscopicMomentum; }

  /** returns the mean potential energy over the macroscopic cell
   * @returns _potentialEnergy potential energy*/
  const double &getPotentialEnergy() const { return _potentialEnergy; }
  /** sets the mean potential energy over the macroscopic cell
   * @param potentialEnergy potential energy*/
  void setPotentialEnergy(const double &potentialEnergy) { _potentialEnergy = potentialEnergy; }

  /** adds a certain amount to the microscopic mass
   * @param mass Mass to be added to the microscopic mass*/
  void addMicroscopicMass(const double &mass) { _microscopicMass += mass; }
  /** adds a certain amount to the microscopic moments
   * @param momentum Momentum to be added to the microscopic moments*/
  void addMicroscopicMomentum(const tarch::la::Vector<dim, double> &momentum) { _microscopicMomentum += momentum; }

  /** adds a certain amount to the macroscopic mass
   * @param mass Mass to be added to the macroscopic mass*/
  void addMacroscopicMass(const double &mass) { _macroscopicMass += mass; }
  /** adds a certain amount to the macroscopic moments
   * @param momentum Momentum to be added to the macroscopic moments*/
  void addMacroscopicMomentum(const tarch::la::Vector<dim, double> &momentum) { _macroscopicMomentum += momentum; }

  /** sets current velocity (sampled right before any distributeX(...) call)
   * @param velocity velocity*/
  void setCurrentVelocity(const tarch::la::Vector<dim, double> &velocity) { _currentVelocity = velocity; }
  /** returns current velocity (sampled right before any distributeX(...) call)
   * @return _currentVelocity velocity*/
  const tarch::la::Vector<dim, double> &getCurrentVelocity() const { return _currentVelocity; }

  /** sets the temperature
   * @param temperature temperature*/
  void setTemperature(const double &temperature) { _temperature = temperature; }
  /** returns the temperature
   * @return _temperature temperature*/
  const double &getTemperature() const { return _temperature; }

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

  /** holds the mean potential energy within the macroscopic cell. This value is
   * needed as a reference
   *  potential energy value for the USHER scheme.
   */
  double _potentialEnergy;

  /** temperature within this cell. Needed for the thermostat. */
  double _temperature;

  /** buffer for current mean velocity in the cell. */
  tarch::la::Vector<dim, double> _currentVelocity;
};

/** describes a quadratic/ cubic macroscopic cell filled with fluid. It is built
 * up by
 *  a certain number of linked cells (from the MD algorithm). The linked
 *  cells need to exactly fill this cell; no overlap/ non-fitting boundaries
 *  shall be allowed.
 *  We can use the MacroscopicCell-structure to evaluate macroscopic quantities
 *  over a certain MD volume and then map macroscopic conserved quantities
 *  between macro- and microscopic simulations.
 *	@brief defines the cell type with cell-averaged quantities. Derived from
 *the class coupling::datastructures::MacroscopicCell
 *	@tparam LinkedCell linked cells that build up the
 * MacroscopicCellWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::datastructures::MacroscopicCellWithLinkedCells : public coupling::datastructures::MacroscopicCell<dim> {
public:
  /** Constructor: initialises the macroscopic cell based on the assumption of
   * having blockSize linked cells;
   *  @param blockSize represents the number of linked cells in all spatial
   * directions.
   */
  MacroscopicCellWithLinkedCells(tarch::la::Vector<dim, unsigned int> blockSize)
      : coupling::datastructures::MacroscopicCell<dim>(), _numberCells(getNumberCells(blockSize)), _linkedCells(NULL) {

    _linkedCells = new LinkedCell *[_numberCells];
    if (_linkedCells == NULL) {
      std::cout << "ERROR coupling::datastructures::MacroscopicCellWithLinkedCells: "
                   "_linkedCells == NULL"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // set each pointer to a NULL pointer
    for (unsigned int i = 0; i < _numberCells; i++) {
      _linkedCells[i] = NULL;
    }
  }
  /** Destructor */
  virtual ~MacroscopicCellWithLinkedCells() {
    if (_linkedCells != NULL) {
      for (unsigned int i = 0; i < _numberCells; i++) {
        _linkedCells[i] = NULL;
      }
      delete[] _linkedCells;
      _linkedCells = NULL;
    }
  }

  /** adds a linked cell to the macroscopic cell and puts it at position index.
We refer to the lexicographic
     *  ordering of the linked cells here.
         * @param cell the linked cell that should be inserted into the
macroscopic cell
         * @param index specifies the position, at which cell shoeld be inserted
     */
  void addLinkedCell(LinkedCell &cell, const unsigned int &index) { _linkedCells[index] = &cell; }

  /** This template fuction applies class A to all linked cells of this
macroscopic cell. The syntax is exactly the same as for regular cell mappings so
     *  that a general cell mapping can also be directly applied to single
macroscopic cells only.
         *	@tparam A
         * 	@param a
     */
  template <class A> void iterateCells(A &a) {
    a.beginCellIteration();
    for (unsigned int i = 0; i < _numberCells; i++) {
      a.handleCell(*(_linkedCells[i]), i);
    }
    a.endCellIteration();
  }

  /** This template function applies class A to all linked cells of this
macroscopic cell. The syntax is exactly the same as for regular cell mappings so
     *  that a general cell mapping can also be directly applied to single
macroscopic cells only.
     *  This method is const, i.e. id does not modify anything except for the
object a. This allows for further optimisations.
         *	@tparam A
         * 	@param a
     */
  template <class A> void iterateConstCells(A &a) const {
    a.beginCellIteration();
    for (unsigned int i = 0; i < _numberCells; i++) {
      a.handleCell(*(_linkedCells[i]), i);
    }
    a.endCellIteration();
  }

private:
  /** computes and returns the number of cells specified by the product of the
entries of blockSize in all dimensions
        *  @param blockSize represents the number of linked cells in all spatial
directions.
        */
  unsigned int getNumberCells(tarch::la::Vector<dim, unsigned int> blockSize) const {
    unsigned int num = 1;
    for (unsigned int d = 0; d < dim; d++) {
      num = num * blockSize[d];
    }
    return num;
  }

  /** total number of linked cells contained in this macroscopic cell */
  const unsigned int _numberCells;

  /** holds pointers to all linked cells that describe the microscopic dynamics
   * within the macroscopic cell */
  LinkedCell **_linkedCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATA_STRUCTURES_MACROSCOPICCELL_H_
