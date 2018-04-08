// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_DELETEPARTICLEMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_DELETEPARTICLEMAPPING_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include "coupling/datastructures/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class DeleteParticleMapping;
  }
}


/** deletes a certain particle from a macroscopic cell.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::DeleteParticleMapping {
  public:
    DeleteParticleMapping(const unsigned int& particle,coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface),
    _particle(particle), _particleCounter(0),
    _deletedMoleculeCopy(tarch::la::Vector<dim,double>(0.0),tarch::la::Vector<dim,double>(0.0),tarch::la::Vector<dim,double>(0.0),0.0){}

    ~DeleteParticleMapping(){}

    void beginCellIteration(){
      _particleCounter = 0;
    }

    void endCellIteration(){}

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      // return, if we already deleted the respective particle
      if (_particleCounter > _particle){
        return;
      }

      // otherwise: loop over particles
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        // if we reached the respective particle, delete it
        if (_particleCounter == _particle){
          const coupling::interface::Molecule<dim> &myMolecule(it->getConst());

          // save copy of particle
          const tarch::la::Vector<dim,double> position = myMolecule.getPosition();
          const tarch::la::Vector<dim,double> velocity = myMolecule.getVelocity();
          const tarch::la::Vector<dim,double> force = myMolecule.getForce();
          const double potentialEnergy = myMolecule.getPotentialEnergy();
          _deletedMoleculeCopy = coupling::datastructures::Molecule<dim>(position,velocity,force,potentialEnergy);

          // delete molecule from MD simulation
          _mdSolverInterface->deleteMoleculeFromMDSimulation(myMolecule, cell);
          _particleCounter++;
          break;
        }
        _particleCounter++;

        it->next();
      }
      delete it;
    }

    coupling::datastructures::Molecule<dim> getDeletedMolecule() const { return _deletedMoleculeCopy; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const unsigned int _particle;
    unsigned int _particleCounter;
    coupling::datastructures::Molecule<dim> _deletedMoleculeCopy;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_DELETEPARTICLEMAPPING_H_
