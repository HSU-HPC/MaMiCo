// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_

#include <iostream>
#include "coupling/interface/Molecule.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeMeanPotentialEnergyMapping;
  }
}


/** computes the mean potential energy over this macroscopic cell.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeMeanPotentialEnergyMapping {
  public:
    ComputeMeanPotentialEnergyMapping(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface), _meanPotentialEnergy(0.0),_particleCounter(0){}

    ~ComputeMeanPotentialEnergyMapping(){}

    void beginCellIteration(){
      _meanPotentialEnergy = 0.0;
      _particleCounter = 0;
    }

    void endCellIteration(){
      if (_particleCounter!=0){
        _meanPotentialEnergy = _meanPotentialEnergy/_particleCounter;
      }
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        const coupling::interface::Molecule<dim> &wrapper(it->getConst());
        _meanPotentialEnergy += wrapper.getPotentialEnergy();
        _particleCounter++;

        it->next();
      }
      delete it;
    }

    double getPotentialEnergy() const { return _meanPotentialEnergy; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    double _meanPotentialEnergy;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMEANPOTENTIALENERGYMAPPING_H_
