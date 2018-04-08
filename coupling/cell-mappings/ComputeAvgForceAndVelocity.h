// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeAvgForceAndVelocity;
  }
}


/** sums up all force/velocity vectors and counts molecules inside a linked cell. Afterwards, the average force/velocity
 *  contribution is computed.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeAvgForceAndVelocity {
  public:
    ComputeAvgForceAndVelocity(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface), _force(0.0),_velocity(0.0),_particleCounter(0){}

    ~ComputeAvgForceAndVelocity(){}

    void beginCellIteration(){
      _particleCounter = 0;
      _force = tarch::la::Vector<dim,double>(0.0);
      _velocity = tarch::la::Vector<dim,double>(0.0);
    }

    void endCellIteration(){
      if (_particleCounter!=0){
        _force = (1.0/_particleCounter)*_force;
        _velocity = (1.0/_particleCounter)*_velocity;
      }
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        coupling::interface::Molecule<dim> &wrapper(it->get());
        _particleCounter++;
        _force = _force + wrapper.getForce();
        _velocity = _velocity + wrapper.getVelocity();
        it->next();
      }
      delete it;
    }

    tarch::la::Vector<dim,double> getAvgForce() const { return _force;}
    tarch::la::Vector<dim,double> getAvgVelocity() const { return _velocity;}

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    tarch::la::Vector<dim,double> _force;
    tarch::la::Vector<dim,double> _velocity;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEAVGFORCEANDVELOCITY_H_
