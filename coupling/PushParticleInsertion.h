// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_PUSHPARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_PUSHPARTICLEINSERTION_H_

#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "coupling/interface/Molecule.h"
#include "coupling/datastructures/DummyCell.h"
#include "coupling/datastructures/Molecule.h"
#include "coupling/ParticleInsertion.h"
#include <iostream>
#include <cstdlib>


namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class PushParticleInsertion;
}

template<class LinkedCell,unsigned int dim>
class Pusher{
public:
  Pusher(coupling::interface::MDSolverInterface<LinkedCell, dim> * const mdSolverInterface,
    const tarch::la::Vector<dim,double>& lowerLeftFront, const tarch::la::Vector<dim,double>& upperRightBack,
    const tarch::la::Vector<dim,double>& toMove, const unsigned int& indexOfDirection, const tarch::la::Vector<dim,double>& meanVelocity,
    const double& temperature): _mdSolverInterface(mdSolverInterface),
     _lowerLeftFront(lowerLeftFront), _upperRightBack(upperRightBack), _toMove(toMove), _indexOfDirection(indexOfDirection),
     _meanVelocity(meanVelocity), _temperature(temperature){}

  void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
    if(cell.hasGhostMolecules()){
      // coupling::interface::MoleculeIterator<LinkedCell,dim> *it = static_cast<coupling::interface::SimpleMDSolverInterface*>(_mdSolverInterface)->getGhostMoleculeIterator(cell);
      coupling::interface::SimpleMDGhostMoleculeIterator *it = static_cast<coupling::interface::SimpleMDSolverInterface*>(_mdSolverInterface)->getGhostMoleculeIterator(cell);
      it->begin(_indexOfDirection);
      while(it->continueIteration(_indexOfDirection)){
        coupling::interface::Molecule<dim> &wrapper(it->get());
        // extract position of each molecule
        tarch::la::Vector<dim,double> position(wrapper.getPosition()+_toMove);
        if(position[_indexOfDirection]>_lowerLeftFront[_indexOfDirection] && position[_indexOfDirection]<_upperRightBack[_indexOfDirection]){
          coupling::datastructures::Molecule<dim> molecule;
          tarch::la::Vector<dim,double> velocity(0.0);
          _mdSolverInterface->getInitialVelocity(_meanVelocity,_mdSolverInterface->getKB(),_temperature,velocity);
          molecule.setVelocity(velocity);
          molecule.setPosition(position);
          _mdSolverInterface->calculateForceAndEnergy(molecule);
          // add molecule to MD simulation and linked cell structures
          _mdSolverInterface->addMoleculeToMDSimulation(molecule);
          position[_indexOfDirection] -= 2.5;
          _insertedParticles++;
        }
        wrapper.setPosition(position);
        it->next();
      }
      delete it;
    }
  }

  void beginCellIteration(){}

  void endCellIteration(){}

  const unsigned int getInsertedParticles(){return _insertedParticles;}

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim> * const _mdSolverInterface;
  const tarch::la::Vector<dim,double> _lowerLeftFront;
  const tarch::la::Vector<dim,double> _upperRightBack;
  const tarch::la::Vector<dim,double> _toMove;
  const unsigned int _indexOfDirection;
  const tarch::la::Vector<dim,double> _meanVelocity;
  const double _temperature;
  unsigned int _insertedParticles{0};
};

template<class LinkedCell,unsigned int dim>
class DeleteOutestParticle {
  public:
    DeleteOutestParticle(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
    const tarch::la::Vector<dim,double> lowerLeftFront, const tarch::la::Vector<dim,double> upperRightBack,
    const unsigned int indexOfDirection): _mdSolverInterface(mdSolverInterface), _lowerLeftFront(lowerLeftFront),
    _upperRightBack(upperRightBack), _indexOfDirection(indexOfDirection){}
    ~DeleteOutestParticle(){}

    void beginCellIteration(){}

    void endCellIteration(){}

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      // std::cout << "4 ";
      it->begin();
      coupling::interface::Molecule<dim>& deleteMolecule(it->get());
      const tarch::la::Vector<dim,double> position = deleteMolecule.getPosition();
      double closestWallPosition;
      double actualPosition = position[_indexOfDirection];
      // std::cout << "5 ";
      if( (_lowerLeftFront[_indexOfDirection]-actualPosition) < (actualPosition-_upperRightBack[_indexOfDirection]) ){
        closestWallPosition = _lowerLeftFront[_indexOfDirection];}
      else{ closestWallPosition = _upperRightBack[_indexOfDirection];}
      it->next();
      // find most outer particle in flow direction
      while(it->continueIteration()){
        // std::cout << "7 ";
        const coupling::interface::Molecule<dim> &myMolecule(it->getConst());
        const tarch::la::Vector<dim,double> position = myMolecule.getPosition();
        if(abs(position[_indexOfDirection]-closestWallPosition) < abs(actualPosition-closestWallPosition) ){
          actualPosition = position[_indexOfDirection];
          deleteMolecule = it->get();
        }
        // std::cout << "8 ";
        it->next();
      }
      // std::cout << "9 ";
      // delete molecule from MD simulation
      _mdSolverInterface->deleteMoleculeFromMDSimulation(deleteMolecule, cell);
      // std::cout << "particle was deleted " << std::endl;
      delete it;
    }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    const tarch::la::Vector<dim,double> _lowerLeftFront;
    const tarch::la::Vector<dim,double> _upperRightBack;
    const unsigned int _indexOfDirection;
};

/** An alternative insertion strategy: Particles will be pushed into the domain.
 *  @author Helene Wittenberg
 */
template<class LinkedCell,unsigned int dim>
class coupling::PushParticleInsertion: public coupling::ParticleInsertion<LinkedCell,dim> {
  public:
    PushParticleInsertion(unsigned int deleteMassEveryTimestep,
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface, const coupling::IndexConversion<dim> &indexConversion):
    ParticleInsertion<LinkedCell, dim>(deleteMassEveryTimestep),
    _mdSolverInterface(mdSolverInterface), _lowerLeftFront(indexConversion.getGlobalMDDomainOffset()),
    _upperRightBack(indexConversion.getGlobalMDDomainSize()+_lowerLeftFront), _moleculeMass(_mdSolverInterface->getMoleculeMass())
    {
      tarch::la::Vector<dim,double> linkedCellSize(2.5, 2.5, 2.5); // get information from somewhere (not indexConversion or mdSolverInterface)

      // left wall
      tarch::la::Vector<dim,double> linkedCellOffset=_lowerLeftFront;
      tarch::la::Vector<dim,double> moveTo(-linkedCellSize[0], 0, 0);
      for(linkedCellOffset[1]=_lowerLeftFront[1]; linkedCellOffset[1]<_upperRightBack[1]; linkedCellOffset[1]+=linkedCellSize[1]){
        for(linkedCellOffset[2]=_lowerLeftFront[2]; linkedCellOffset[2]<_upperRightBack[2]; linkedCellOffset[2]+=linkedCellSize[2]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo,1);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo,4);
        }
      }
      // front wall
      linkedCellOffset=_lowerLeftFront;
      moveTo = tarch::la::Vector<3,double>(0, -linkedCellSize[1], 0);
      for(linkedCellOffset[0]=_lowerLeftFront[0]; linkedCellOffset[0]<_upperRightBack[0]; linkedCellOffset[0]+=linkedCellSize[0]){
        for(linkedCellOffset[2]=_lowerLeftFront[2]; linkedCellOffset[2]<_upperRightBack[2]; linkedCellOffset[2]+=linkedCellSize[2]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo, 2);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo, 5);
        }
      }
      // lower wall
      linkedCellOffset=_lowerLeftFront;
      moveTo = tarch::la::Vector<3,double>(0, 0, -linkedCellSize[2]);
      for(linkedCellOffset[1]=_lowerLeftFront[1]; linkedCellOffset[1]<_upperRightBack[1]; linkedCellOffset[1]+=linkedCellSize[1]){
        for(linkedCellOffset[0]=_lowerLeftFront[0]; linkedCellOffset[0]<_upperRightBack[0]; linkedCellOffset[0]+=linkedCellSize[0]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo, 3);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo, 6);
        }
      }
      // right wall
      linkedCellOffset=_upperRightBack-linkedCellSize;
      moveTo = tarch::la::Vector<dim,double>(linkedCellSize[0], 0, 0);
      for(linkedCellOffset[1]=_upperRightBack[1]-linkedCellSize[1]; linkedCellOffset[1]>=_lowerLeftFront[1]; linkedCellOffset[1]-=linkedCellSize[1]){
        for(linkedCellOffset[2]=_upperRightBack[2]-linkedCellSize[2]; linkedCellOffset[2]>=_lowerLeftFront[2]; linkedCellOffset[2]-=linkedCellSize[2]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo, 1);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo,4);
        }
      }
      // back wall
      linkedCellOffset=_upperRightBack-linkedCellSize;
      moveTo = tarch::la::Vector<3,double>(0, linkedCellSize[1], 0);
      for(linkedCellOffset[0]=_upperRightBack[0]-linkedCellSize[0]; linkedCellOffset[0]>=_lowerLeftFront[0]; linkedCellOffset[0]-=linkedCellSize[0]){
        for(linkedCellOffset[2]=_upperRightBack[2]-linkedCellSize[2]; linkedCellOffset[2]>=_lowerLeftFront[2]; linkedCellOffset[2]-=linkedCellSize[2]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo, 2);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo,5);
        }
      }
      // upper wall
      linkedCellOffset=_upperRightBack-linkedCellSize;
      moveTo = tarch::la::Vector<3,double>(0, 0, linkedCellSize[2]);
      for(linkedCellOffset[1]=_upperRightBack[1]-linkedCellSize[1]; linkedCellOffset[1]>=_lowerLeftFront[1]; linkedCellOffset[1]-=linkedCellSize[1]){
        for(linkedCellOffset[0]=_upperRightBack[0]-linkedCellSize[0]; linkedCellOffset[0]>=_lowerLeftFront[0]; linkedCellOffset[0]-=linkedCellSize[0]){
          coupling::datastructures::DummyCell newMolecules(linkedCellOffset);
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[0], moveTo, 3);
          for(unsigned int i=1; i<newMolecules.size()-1; i++){
            static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[i], moveTo);
          }
          static_cast<coupling::interface::SimpleMDSolverInterface*>(mdSolverInterface)->addGhostMoleculeToMDSimulation(newMolecules[newMolecules.size()-1], moveTo,6);
        }
      }
    }

    ~PushParticleInsertion()override{}

    typename coupling::ParticleInsertion<LinkedCell,dim>::Action insertDeleteMass(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const tarch::la::Vector<dim,double>& macroscopicCellPosition,
      const tarch::la::Vector<dim,double>& macroscopicCellSize,
      const tarch::la::Vector<dim,double>& meanVelocity, const double &temperature,
      const coupling::BoundaryForceController<LinkedCell,dim>& boundaryForceController
    ) override{
      typename coupling::ParticleInsertion<LinkedCell, dim>::Action resultAction = coupling::ParticleInsertion<LinkedCell, dim>::NoAction;

      if(cell.getMassFlux()>0.0){
        const unsigned int insertions = moveGhostParticles(cell, meanVelocity, temperature);
        if(insertions>0){
          resultAction=coupling::ParticleInsertion<LinkedCell, dim>::Insertion;
          cell.addMicroscopicMass(-static_cast<float>(insertions)*_moleculeMass);
        }
      }
      else if(cell.getMicroscopicMass()<-0.6*_moleculeMass){
        resultAction = deleteMass(cell);
        if(resultAction==coupling::ParticleInsertion<LinkedCell, dim>::Deletion){
          cell.addMicroscopicMass(_moleculeMass);
        }
      }
      return resultAction;
    }

    const unsigned int moveGhostParticles(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell, const tarch::la::Vector<dim,double>& meanVelocity, const double &temperature){
      // calculate the push length
      const tarch::la::Vector<dim,double> toMove = cell.getMassFluxVector()*_densityFactor;
      // get the direction (dimension) for the mass insertion
      unsigned int indexOfDirection = 0;
      while(toMove[indexOfDirection]==0)
        indexOfDirection++;
      // push the particles into the domain
      Pusher<LinkedCell,dim> pusher(_mdSolverInterface, _lowerLeftFront, _upperRightBack, toMove, indexOfDirection, meanVelocity, temperature);
      cell.iterateCells(pusher);
      return pusher.getInsertedParticles();
    }

    typename coupling::ParticleInsertion<LinkedCell,dim>::Action deleteMass(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell){
      // count particles with computeMassMapping
      coupling::cellmappings::ComputeMassMapping<LinkedCell,dim> computeMassMapping(_mdSolverInterface);
      cell.iterateConstCells(computeMassMapping);
      if(computeMassMapping.getNumberOfParticles()<1){
        return coupling::ParticleInsertion<LinkedCell,dim>::NoAction;
      }
      const tarch::la::Vector<dim,double> massFlux = cell.getMassFluxVector();
      unsigned int indexOfDirection = 0;
      while(massFlux[indexOfDirection]==0)
        indexOfDirection++;
      DeleteOutestParticle<LinkedCell,dim> deleteOutestParticle(_mdSolverInterface, _lowerLeftFront, _upperRightBack, indexOfDirection);
      cell.iterateCells(deleteOutestParticle);
      return coupling::ParticleInsertion<LinkedCell,dim>::Deletion;
    }

    /** The push algortihm doesn't require the potential energy landscape, therefore false is returned */
    bool requiresPotentialEnergyLandscape()override{return false;}

  private:
    coupling::interface::MDSolverInterface<LinkedCell, dim> * const _mdSolverInterface;
    const tarch::la::Vector<dim,double> _lowerLeftFront;
    const tarch::la::Vector<dim,double> _upperRightBack;
    const double _moleculeMass;
    const double _densityFactor{1/2.5/2.5/0.768}; // density of the cells with the ghost molecules

};
#endif // _MOLECULARDYNAMICS_COUPLING_PUSHPARTICLEINSERTION_H_
