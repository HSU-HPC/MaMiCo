// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_

#include <iostream>
#include <list>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeKineticEnergyMapping;
  }
}


/** computes the kinetic energy.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeKineticEnergyMapping {
  public:
    ComputeKineticEnergyMapping(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): _mdSolverInterface(mdSolverInterface), _kineticEnergy(0.0){}

    ~ComputeKineticEnergyMapping(){}

    void beginCellIteration(){
      _kineticEnergy = 0.0;
    }

    void endCellIteration(){
      _kineticEnergy = 0.5*_mdSolverInterface->getMoleculeMass()*_kineticEnergy;
    }

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){

      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        const coupling::interface::Molecule<dim> &wrapper(it->getConst());
        _kineticEnergy += tarch::la::dot(wrapper.getVelocity(),wrapper.getVelocity());

        it->next();
      }
      delete it;
    }

    double getKineticEnergy() const { return _kineticEnergy; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    double _kineticEnergy;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEKINETICENERGYMAPPING_H_
