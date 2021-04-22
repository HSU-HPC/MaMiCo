// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_LIMITFORCE_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_LIMITFORCE_H_

#include <cmath>
#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class LimitForce;
  }
}

/**
 *
 *  @author Helene Wittenberg
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::LimitForce {
  public:
    LimitForce(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):_mdSolverInterface(mdSolverInterface){}

    ~LimitForce(){}

    void beginCellIteration(){}

    void endCellIteration(){}

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        coupling::interface::Molecule<dim> &wrapper(it->get());

        // extract position and forc of each molecule
        tarch::la::Vector<dim,double> force(wrapper.getForce());
        for(unsigned int i = 0; i<dim; i++){
          if(force[i]>_maxForce){
            double factor = _maxForce/force[i];
            force[0]=force[0]*factor;
            force[1]=force[1]*factor;
            force[2]=force[2]*factor;
          }
          else if(force[i]<-_maxForce){
            double factor = -_maxForce/force[i];
            force[0]=force[0]*factor;
            force[1]=force[1]*factor;
            force[2]=force[2]*factor;
          }
        }
        wrapper.setForce(force);
        it->next();
      }
      delete it;
    }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const double _maxForce{500.0};
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_LIMITFORCE_H_
