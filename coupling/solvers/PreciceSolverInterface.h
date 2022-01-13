// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVERINTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
  namespace solvers {
    template<unsigned int dim>
    class PreciceSolverInterface;
  }
}


/** 
 *  @author Louis Viot
 */
template<unsigned int dim>
class coupling::solvers::PreciceSolverInterface: public coupling::interface::MacroscopicSolverInterface<dim> {
public:
PreciceSolverInterface(tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,unsigned int outerRegion=1) :
coupling::interface::MacroscopicSolverInterface<dim>(),
_outerRegion(outerRegion),
_globalNumberMacroscopicCells(globalNumberMacroscopicCells)
{
}

~PreciceSolverInterface()
{
}

bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim,unsigned int> globalCellIndex) override 
{
	bool recv=true;
  	for (unsigned int d = 0; d < dim; d++)
 	{ 
  		recv &= (globalCellIndex[d]>_outerRegion-1) && (globalCellIndex[d]<_globalNumberMacroscopicCells[d]+2-_outerRegion); 
	}
	return recv;
}


bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim,unsigned int> globalCellIndex) override
{
	bool outer=false;
	bool recv=true;
  	for (unsigned int d = 0; d < dim; d++)
  	{ 
  		recv &= (globalCellIndex[d]>_outerRegion) && (globalCellIndex[d]<_globalNumberMacroscopicCells[d]+1-_outerRegion); 
	}
  	for (unsigned int d = 0; d < dim; d++)
  	{ 
  		outer|= (globalCellIndex[d]<1) || (globalCellIndex[d]>_globalNumberMacroscopicCells[d]); 
	}
	return (!outer) && (!recv);
}

std::vector<unsigned int> getRanks(tarch::la::Vector<dim,unsigned int> globalCellIndex)
{
	std::vector<unsigned int> ranks; 
	ranks.push_back(0); 
	return ranks; 
}

private:
	const unsigned int _outerRegion;
	const tarch::la::Vector<dim,unsigned int> _globalNumberMacroscopicCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVERINTERFACE_H_
