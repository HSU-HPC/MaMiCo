// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELL_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_MACROSCOPICCELL_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "tarch/la/Vector.h"

namespace coupling {
  namespace datastructures {

    template<unsigned int dim>
    class MacroscopicCell;
    template<class LinkedCell,unsigned int dim>
    class MacroscopicCellWithLinkedCells;
  }
}



/** cell type with cell-averaged quantities only (no linked cells).
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::datastructures::MacroscopicCell {
  public:
    /** initialises the macroscopic cell based on the assumption of having blockSize linked cells;
     *  blockSize represents the number of linked cells in all spatial directions.
     */
    MacroscopicCell():
      _microscopicMass(0.0),_microscopicMomentum(0.0),
      _macroscopicMass(0.0),_macroscopicMomentum(0.0),
      _potentialEnergy(0.0),_temperature(0.0),_currentVelocity(0.0){}

    virtual ~MacroscopicCell(){}

    /** sets/ returns the microscopic moments */
    void setMicroscopicMass(const double& mass){ _microscopicMass = mass; }
    const double& getMicroscopicMass() const {return _microscopicMass;}
    void setMicroscopicMomentum(const tarch::la::Vector<dim,double>& momentum){ _microscopicMomentum = momentum; }
    const tarch::la::Vector<dim,double>& getMicroscopicMomentum() const {return _microscopicMomentum;}

    /** sets/ returns the macroscopic moments */
    void setMacroscopicMass(const double& mass){ _macroscopicMass = mass; }
    const double& getMacroscopicMass() const {return _macroscopicMass;}
    void setMacroscopicMomentum(const tarch::la::Vector<dim,double>& momentum){ _macroscopicMomentum = momentum; }
    const tarch::la::Vector<dim,double>& getMacroscopicMomentum() const {return _macroscopicMomentum;}

    /** sets/ returns the mean potential energy over the macroscopic cell */
    const double& getPotentialEnergy() const { return _potentialEnergy; }
    void setPotentialEnergy(const double& potentialEnergy){ _potentialEnergy=potentialEnergy;}


    /** adds a certain amount to the microscopic moments*/
    void addMicroscopicMass(const double& mass){ _microscopicMass += mass;}
    void addMicroscopicMomentum(const tarch::la::Vector<dim,double>& momentum){ _microscopicMomentum += momentum; }

    /** adds a certain amount to the macroscopic moments*/
    void addMacroscopicMass(const double& mass){ _macroscopicMass += mass;}
    void addMacroscopicMomentum(const tarch::la::Vector<dim,double>& momentum){ _macroscopicMomentum += momentum; }

    /** sets/ gets current velocity (sampled right before any distributeX(...) call) */
    void setCurrentVelocity(const tarch::la::Vector<dim,double>& velocity){ _currentVelocity = velocity;}
    const tarch::la::Vector<dim,double>& getCurrentVelocity() const {return _currentVelocity;}

    /** sets and return the temperature value */
    void setTemperature(const double& temperature){ _temperature = temperature;}
    const double& getTemperature() const { return _temperature;}

    void setMassFlux(const double& massFlux){_massFlux = massFlux;}
    const double& getMassFlux() const {return _massFlux;}


    /** buffers for macroscopic quantities of mass, momentum that need to be transferred
     *  from the macroscopic to the microscopic simulation */
    double _microscopicMass;
    tarch::la::Vector<dim,double> _microscopicMomentum;

    /** buffers for macroscopic quantities of mass, momentum that need to be returned
     *  from the microscopic to the macroscopic simulation */
    double _macroscopicMass;
    tarch::la::Vector<dim,double> _macroscopicMomentum;

    /** holds the mean potential energy within the macroscopic cell. This value is needed as a reference
     *  potential energy value for the USHER scheme.
     */
    double _potentialEnergy;

    /** temperature within this cell. Needed for the thermostat. */
    double _temperature;

    /** buffer for current mean velocity in the cell. */
    tarch::la::Vector<dim,double> _currentVelocity;

    double _massFlux{0.0};
};


/** describes a quadratic/ cubic macroscopic cell filled with fluid. It is built up by
 *  a certain number of linked cells (from the MD algorithm). The linked
 *  cells need to exactly fill this cell; no overlap/ non-fitting boundaries
 *  shall be allowed.
 *  We can use the MacroscopicCell-structure to evaluate macroscopic quantities
 *  over a certain MD volume and then map macroscopic conserved quantities
 *  between macro- and microscopic simulations.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::datastructures::MacroscopicCellWithLinkedCells: public coupling::datastructures::MacroscopicCell<dim> {
  public:
    /** initialises the macroscopic cell based on the assumption of having blockSize linked cells;
     *  blockSize represents the number of linked cells in all spatial directions.
     */
    MacroscopicCellWithLinkedCells(tarch::la::Vector<dim,unsigned int> blockSize):
      coupling::datastructures::MacroscopicCell<dim>(),
      _numberCells(getNumberCells(blockSize)),
      _linkedCells(NULL){

      _linkedCells = new LinkedCell* [_numberCells];
      if (_linkedCells == NULL){
        std::cout << "ERROR coupling::datastructures::MacroscopicCellWithLinkedCells: _linkedCells == NULL" << std::endl;
        exit(EXIT_FAILURE);
      }
      // set each pointer to a NULL pointer
      for (unsigned int i = 0; i < _numberCells; i++){
        _linkedCells[i] = NULL;
      }
    }
    virtual ~MacroscopicCellWithLinkedCells(){
      if (_linkedCells != NULL){
        for (unsigned int i = 0; i < _numberCells; i++){ _linkedCells[i] = NULL;}
        delete [] _linkedCells;
        _linkedCells = NULL;
      }
    }

    /** adds a linked cell to the macroscopic cell and puts it at position index. We refer to the lexicographic
     *  ordering of the linked cells here.
     */
    void addLinkedCell(LinkedCell& cell, const unsigned int& index){
      _linkedCells[index] = &cell;
    }

    /** applies class A to all linked cells of this macroscopic cell. The syntax is exactly the same as for regular cell mappings so
     *  that a general cell mapping can also be directly applied to single macroscopic cells only.
     */
    template<class A>
    void iterateCells(A& a){
      a.beginCellIteration();
      for (unsigned int i = 0; i < _numberCells; i++){
         a.handleCell( *(_linkedCells[i]),i);
      }
      a.endCellIteration();
    }

    /** applies class A to all linked cells of this macroscopic cell. The syntax is exactly the same as for regular cell mappings so
     *  that a general cell mapping can also be directly applied to single macroscopic cells only.
     *  This method is const, i.e. id does not modify anything except for the object a. This allows for further
     *  optimisations.
     */
    template<class A>
    void iterateConstCells(A& a) const{
      a.beginCellIteration();
      for (unsigned int i = 0; i < _numberCells; i++){
         a.handleCell( *(_linkedCells[i]),i);
      }
      a.endCellIteration();
    }

  private:
    /** computes and returns the number of cells specified by the product of the entries of blockSize */
    unsigned int getNumberCells(tarch::la::Vector<dim,unsigned int> blockSize) const {
      unsigned int num = 1;
      for (unsigned int d = 0; d < dim; d++){
        num = num*blockSize[d];
      }
      return num;
    }

    /** total number of linked cells contained in this macroscopic cell */
    const unsigned int _numberCells;

    /** holds pointers to all linked cells that describe the microscopic dynamics within the macroscopic cell */
    LinkedCell** _linkedCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_DATA_STRUCTURES_MACROSCOPICCELL_H_
