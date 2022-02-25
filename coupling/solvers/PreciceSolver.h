#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "coupling/solvers/CouetteSolver.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/IndexConversion.h"
#include "precice/SolverInterface.hpp"

namespace coupling {
  namespace solvers{
    class PreciceSolver;
  }
}

class coupling::solvers::PreciceSolver: public coupling::solvers::AbstractCouetteSolver<3> {
public:
	PreciceSolver(
	  const double channelHeight,
	  const double dx,
	  const double dt,
	  const int plotEveryTimestep,
	  const std::string filestem) :
	coupling::solvers::AbstractCouetteSolver<3>(),
	_channelHeight(channelHeight),
	_dx(dx), 
	_dt(dt), 
	_plotEveryTimestep(plotEveryTimestep), 
	_filestem(filestem)
	{
		std::cout << "PreciceSolver constructor call" << std::endl;	
		_interface = new precice::SolverInterface("MAMICO","precice-config.xml",0,1);
		int dim = _interface->getDimensions();
		if (dim!=3)
		{
			std::cout << "main: preCICE dimension should be 3" << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "_channelHeight " << _channelHeight << std::endl;
		std::cout << "_dx " << _dx << std::endl;
		std::cout << "_domainSize " << _domainSize << std::endl;		
	}

    virtual ~PreciceSolver()
    {
		if (_vel !=NULL){delete [] _vel; _vel=NULL;}
		if (_density!=NULL){delete [] _density; _density=NULL;}
		if (_interface!=NULL){_interface->finalize(); delete _interface; _interface=NULL;}
    }
    
  	void init(	tarch::la::Vector<3, double> mdDomainOffset, 
  				tarch::la::Vector<3, double> mdDomainSize, 
  				unsigned int overlapStrip,
  				const coupling::IndexConversion<3> &indexConversion) 
  	{
	  	const tarch::la::Vector<3, unsigned int> cells = indexConversion.getGlobalNumberMacroscopicCells();
		const unsigned int numberOfCells = cells[0] * cells[1] * cells[2];
		for (unsigned int i = 0; i < numberOfCells; i++) {
			const tarch::la::Vector<3, unsigned int> cellVectorIndex = indexConversion.getGlobalVectorCellIndex(i);
			tarch::la::Vector<3, double> coords;
			for (int d = 0; d < 3; d++) {
				coords[d]=mdDomainOffset[d] + cellVectorIndex[d]*_dx;
			}		
		}
		int meshID = _interface->getMeshID("MamicoMesh");
		double* coords = new double[dim*_numberOfVertices];
		_vel = new double[dim*_numberOfVertices];
		_density = new double[_numberOfVertices];
		for (int vertexIndex = 0; vertexIndex < _numberOfVertices; vertexIndex++)
		{
			for (int currentDim = 0; currentDim < dim; currentDim++)
				_vel[dim*vertexIndex+currentDim]  = (double)0.0;
			_density[vertexIndex] = 1.0;
		}
		_vertexIDs = new int[_numberOfVertices];
		_interface->setMeshVertices(meshID, _numberOfVertices, coords, _vertexIDs); 
		delete[] coords;
		_precice_dt = _interface->initialize();
		std::cout << "PreciceSolver constructor called" << std::endl;
		int exvsID = _interface->getDataID("ExternalVelocities", meshID); 
		double* velocities = new double[_numberOfVertices*dim];
		_interface->readBlockVectorData(exvsID, _numberOfVertices, _vertexIDs, velocities);
		for (int vertexIndex = 0; vertexIndex < _numberOfVertices; vertexIndex++)
		{
			std::cout << "PreciceSolver: read velocities[" << vertexIndex << "] = [" << velocities[dim*vertexIndex] << "," << velocities[dim*vertexIndex+1] << "," << velocities[dim*vertexIndex+2] << "]" << std::endl;
		}
	}
    
	void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override
	{
		_wallVelocity = wallVelocity;
	}
	
	tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
	{
		tarch::la::Vector<3,unsigned int> coords;
		for (unsigned int d = 0; d < 3; d++) 
			coords[d] =  (unsigned int) pos[d]/_dx;
		const int index = get(coords[0],coords[1],coords[2]);
		tarch::la::Vector<3,double> vel(0.0);
		for (int d = 0; d < 3; d++)
			vel[d] = _vel[3*index+d];
		return vel;
	}
	
	double getDensity(tarch::la::Vector<3,double> pos) const
	{
		tarch::la::Vector<3,unsigned int> coords;
		for (unsigned int d = 0; d < 3; d++)
			coords[d] =  (unsigned int) pos[d]/_dx;
		const int index = get(coords[0],coords[1],coords[2]);
		return _density[index];
	}
	
	void advance(double dt) override 
	{
		if(skipRank()){return;}
		std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
		if (_interface->isCouplingOngoing()) {
			int dim = _interface->getDimensions();
			int meshID = _interface->getMeshID("MamicoMesh");
			int exvsID = _interface->getDataID("ExternalVelocities", meshID); 
			double* velocities = new double[_numberOfVertices*dim];
			_interface->readBlockVectorData(exvsID, _numberOfVertices, _vertexIDs, velocities);
			for (int vertexIndex = 0; vertexIndex < _numberOfVertices; vertexIndex++)
			{
				std::cout << "PreciceSolver::advance: read velocities[" << vertexIndex << "] = [" << velocities[dim*vertexIndex] << "," << velocities[dim*vertexIndex+1] << "," << velocities[dim*vertexIndex+2] << "]" << std::endl;
			}
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
			std::cout << "PreciceSolver::advance preCICE solver advanced with dt = " << computed_dt << std::endl;
    	}
	}

private:    
    const double _channelHeight; //
    const double _dx;
    const double _dt;
    const int _plotEveryTimestep;
    const std::string _filestem;
    const int _domainSize;
    const int _numberOfVertices;
    int _counter{0};
    double *_vel{NULL};
    double *_density{NULL};
	tarch::la::Vector<3,double> _wallVelocity;
    precice::SolverInterface *_interface;
	double _precice_dt;
	int *_vertexIDs;
	double *_vertexCoords;
};

#endif
