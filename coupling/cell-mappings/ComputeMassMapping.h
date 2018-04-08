// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMASSMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMASSMAPPING_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeMassMapping;
  }
}


/** computes the mass over certain linked cells.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeMassMapping {
  public:
    ComputeMassMapping(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface), _mass(0.0),_particleCounter(0){}

    ~ComputeMassMapping(){}

    void beginCellIteration(){
      _mass = 0.0;
      _particleCounter = 0;
    }

    void endCellIteration(){
      _mass = _mdSolverInterface->getMoleculeMass()*_particleCounter;
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        _particleCounter++;
        it->next();
      }
      delete it;
    }

    double getMass() const { return _mass; }
    unsigned int getNumberOfParticles() const { return _particleCounter; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    double _mass;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMASSMAPPING_H_
