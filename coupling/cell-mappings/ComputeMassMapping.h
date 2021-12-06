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


/** 
 *	@brief This class computes the mass over certain linked cells.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::ComputeMassMapping {
  public:
    /** Constructor
	 *	@param mdSolverInterface
	 */
	ComputeMassMapping(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface), _mass(0.0),_particleCounter(0){}

	/** Destructor */
	~ComputeMassMapping(){}

    /** sets the mass and the particlee counter to zero, before the iteration process begins.
	 */
	void beginCellIteration(){
      _mass = 0.0;
      _particleCounter = 0;
    }

    /** computes the mass in a linked cell, by multiplying the number of particles inside the cell with the particel mass(which is assuemd to be constant).
	 */
	void endCellIteration(){
      _mass = _mdSolverInterface->getMoleculeMass()*_particleCounter;
    }

    /** counts the molecules inside a linked cell.
	 *	@param cell
	 *	@param cellIndex
	 */	
	void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        _particleCounter++;
        it->next();
      }
      delete it;
    }

    /** returns the mass inside a linked cell
	 *	@return _mass
	 */
	double getMass() const { return _mass; }
	/** returns the number of particles inside a linked cell
	 *	@return _particleCounter
	 */
    unsigned int getNumberOfParticles() const { return _particleCounter; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    double _mass;
    unsigned int _particleCounter;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_COMPUTEMASSMAPPING_H_
