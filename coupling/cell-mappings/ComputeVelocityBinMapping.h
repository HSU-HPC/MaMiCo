// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEVELOCITYBINMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEVELOCITYBINMAPPING_H_

#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class ComputeVelocityBinMapping;
  }
}

/**
 *	@brief This class computes the momentum over certain linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Helene Wittenberg
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeVelocityBinMapping {
  public:
    /** Constructor
	 *	@param mdSolverInterface
	 */
	ComputeVelocityBinMapping(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
    tarch::la::Vector<dim, double>* velocity, tarch::la::Vector<dim, double>* velocityPlus, tarch::la::Vector<dim, double>* velocityMinus,
    unsigned int* numberParticles, unsigned int* numberParticlesPlus, unsigned int* numberParticlesMinus):
    _mdSolverInterface(mdSolverInterface),_velocity(velocity), _velocityPlus(velocityPlus), _velocityMinus(velocityMinus),
    _numberParticles(numberParticles), _numberParticlesPlus(numberParticlesPlus), _numberParticlesMinus(numberParticlesMinus){}

    /** Destructor */
	~ComputeVelocityBinMapping(){}

    /** sets the mean velocity, momentum and the particle counter to zero, before the iteration process begins.
	 */
	void beginCellIteration(){
    }

    /** computes the mean velocity, momentum in a linked cell,
	 *	by dividing and multiplying the summation of the velocities computed in handleCell(...) over the number if particle and in particle mass respectively.
	 */
	void endCellIteration(){
    }

    /** counts the molecules inside a linked cell and sums up the of the velocity of all particles inside the cell and saves it as momentum.
	 *	@param cell
	 *	@param cellIndex
	 */
	void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
    coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while(it->continueIteration()){
      const coupling::interface::Molecule<dim> &wrapper(it->getConst());
      const int posZ = std::floor(wrapper.getPosition()[2]*4);
      _velocity[posZ] += wrapper.getVelocity();
      _numberParticles[posZ]++;
      if(wrapper.getVelocity()[2]>0.0){
        _velocityPlus[posZ] += wrapper.getVelocity();
        _numberParticlesPlus[posZ]++;
      }
      else{
        _velocityMinus[posZ] += wrapper.getVelocity();
        _numberParticlesMinus[posZ]++;
      }
      it->next();
    }
    delete it;
  }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    tarch::la::Vector<dim, double>* _velocity;
    tarch::la::Vector<dim, double>* _velocityPlus;
    tarch::la::Vector<dim, double>* _velocityMinus;
    unsigned int* _numberParticles;
    unsigned int* _numberParticlesPlus;
    unsigned int* _numberParticlesMinus;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEVELOCITYBINMAPPING_H_
