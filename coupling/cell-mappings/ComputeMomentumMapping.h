// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeMomentumMapping;
  }
}


/** computes the momentum over certain linked cells.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeMomentumMapping {
  public:
    ComputeMomentumMapping(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface),_momentum(0.0),_meanVelocity(0.0),_particleCounter(0){}

    ~ComputeMomentumMapping(){}

    void beginCellIteration(){
      _momentum = tarch::la::Vector<dim,double>(0.0);
      _meanVelocity = tarch::la::Vector<dim,double>(0.0);
      _particleCounter = 0;
    }

    void endCellIteration(){
      if (_particleCounter != 0){
        _meanVelocity = (1.0/((double) _particleCounter))*_momentum;
        _momentum = _mdSolverInterface->getMoleculeMass()*_momentum;
      }
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        const coupling::interface::Molecule<dim> &wrapper(it->getConst());
        _momentum += wrapper.getVelocity();
        _particleCounter++;

        it->next();
      }
      delete it;
    }

    tarch::la::Vector<dim,double> getMomentum() const { return _momentum; }

    tarch::la::Vector<dim,double> getMeanVelocity() const { return _meanVelocity; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    tarch::la::Vector<dim,double> _momentum;
    tarch::la::Vector<dim,double> _meanVelocity;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMOMENTUMMAPPING_H_
