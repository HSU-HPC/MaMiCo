// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VTKMOLECULEPLOTTER_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VTKMOLECULEPLOTTER_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class VTKMoleculePlotter;
  }
}


/** writes molecule data to streams.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::VTKMoleculePlotter {
  public:
    VTKMoleculePlotter(
      std::stringstream &moleculeVelocities,
      std::stringstream &moleculePositions,
      std::stringstream &moleculePotentials,
      const std::string &appendFloatZeros,
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): _mdSolverInterface(mdSolverInterface),
    _moleculeVelocities(moleculeVelocities), _moleculePositions(moleculePositions),
    _moleculePotentials(moleculePotentials), _appendFloatZeros(appendFloatZeros),
    _particleCounter(0){}

    ~VTKMoleculePlotter(){}

    void beginCellIteration(){
      _particleCounter = 0;
    }

    void endCellIteration(){}

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        const coupling::interface::Molecule<dim> &wrapper(it->getConst());
        const tarch::la::Vector<dim,double> position = wrapper.getPosition();
        const tarch::la::Vector<dim,double> velocity = wrapper.getVelocity();
        //std::cout << "Touch molecule " << position << std::endl;
        for (unsigned int d = 0; d < dim; d++){
          _moleculePositions << position[d] << " ";
          _moleculeVelocities << velocity[d] << " ";
        }
        _moleculePositions << _appendFloatZeros << std::endl;
        _moleculeVelocities << _appendFloatZeros << std::endl;
        _moleculePotentials << wrapper.getPotentialEnergy() << std::endl;

        _particleCounter++;
        it->next();
      }
      delete it;
    }

    const unsigned int& getParticleCounter() const { return _particleCounter;}

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    std::stringstream &_moleculeVelocities;
    std::stringstream &_moleculePositions;
    std::stringstream &_moleculePotentials;
    const std::string &_appendFloatZeros;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_VTKMOLECULEPLOTTER_H_
