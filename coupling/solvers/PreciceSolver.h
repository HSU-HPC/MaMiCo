// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_SOLVERS_PRECICESOLVER_H_
#define _COUPLING_SOLVERS_PRECICESOLVER_H_

#include "coupling/solvers/CouetteSolver.h"
#include "precice/SolverInterface.hpp"

namespace coupling{
  namespace solvers{
    class PreciceSolver;
  }
}

/**
 * @author Louis Viot
 */
class coupling::solvers::PreciceSolver: public coupling::solvers::AbstractCouetteSolver<3> {
public:
PreciceSolver(int rank, int plotEveryTimestep, double channelheight): 
AbstractCouetteSolver<3>(),
_channelHeight(channelheight),
_rank(rank),
_plotEveryTimestep(plotEveryTimestep),
_timeStepCounter(0)
{
	if(skipRank()){return;}
}

virtual ~PreciceSolver()
{
	if(skipRank()){return;}
}


void advance(double dt) override
{
	if(skipRank()){return;}
    std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
	plottxt();	
	_timeStepCounter++;
}

tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
{
	return tarch::la::Vector<3,double>(0, 0, 0);
};

void setWallVelocity(tarch::la::Vector<3,double> wallVelocity) override
{
};
  	
const tarch::la::Vector<3,double> getOuterPointFromBoundary(const int layer, const int index)
{
	return tarch::la::Vector<3,double>(0, 0, 0);
}

private:
void plottxt() 
{
	if (_plotEveryTimestep < 1 || _timeStepCounter % _plotEveryTimestep > 0) return;
	std::stringstream ss; 
	ss << "velocity_" << _timeStepCounter << ".txt";
	std::ofstream file(ss.str().c_str());
	if (!file.is_open())
	{
		std::cout << "ERROR PreciceSolver::plottxt(): Could not open file " << ss.str() << "!" << std::endl; 
		exit(EXIT_FAILURE);
	}
	std::stringstream velocity;

	// loop over domain (incl. boundary)
	double _dx=_channelHeight/10;
	double y=_channelHeight/2;
	double x=_channelHeight/2;
	for (double z = 0.0; z < _channelHeight; z+=_dx)
	{
		velocity << z << ", " << z << ", " << z << std::endl;
	}
	file << velocity.str() << std::endl;
	file.close();
}

bool skipRank(){
	return _rank!=0;
}

double _channelHeight;
int _rank;
int _plotEveryTimestep;
int _timeStepCounter;
};
#endif // _COUPLING_SOLVERS_PRECICESOLVER_H
