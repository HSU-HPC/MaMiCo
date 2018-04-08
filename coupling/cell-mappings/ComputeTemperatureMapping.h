// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeTemperatureMapping;
  }
}


/** computes the temperature in a certain (macroscopic) cell.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeTemperatureMapping {
  public:
    ComputeTemperatureMapping(const tarch::la::Vector<dim,double>& meanVelocity,coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface),
    _meanVelocity(meanVelocity),
    _temperature(0.0), _particleCounter(0){}

    ~ComputeTemperatureMapping(){}

    void beginCellIteration(){
      _temperature = 0.0;
      _particleCounter = 0;
    }

    void endCellIteration(){
      _temperature = _temperature * _mdSolverInterface->getMoleculeMass();
      if (_particleCounter != 0){
        _temperature = _temperature/ (dim*_mdSolverInterface->getKB()*_particleCounter);
      }
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        const coupling::interface::Molecule<dim> &wrapper(it->getConst());
        _temperature += tarch::la::dot((wrapper.getVelocity() - _meanVelocity),(wrapper.getVelocity() - _meanVelocity));
        _particleCounter++;

        it->next();
      }
      delete it;
    }

    double getTemperature() const { return _temperature; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const tarch::la::Vector<dim,double> _meanVelocity;
    double _temperature;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTETEMPERATUREMAPPING_H_
