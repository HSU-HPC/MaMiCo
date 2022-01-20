// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
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
	/** @brief a simple constructor
	 *  @param rank of this solver 
	 *  @param channelheight the width and height of the channel in y and z direction
	 *  @param dx the spacial step size, and equidistant grid is applied
	 *  @param dt the time step
	 *  @param kinVisc the kinematic viscosity of the fluid
	 *  @param plotEveryTimestep the time step interval for plotting data;
	 *                           4 means, every 4th time step is plotted
	 *  @param filestem the name of the plotted file
	 */
	PreciceSolver(
	  const int rank,
	  const double channelHeight,
	  const double dx,
	  const double dt,
	  const double kinVisc,
	  const int plotEveryTimestep,
	  const std::string filestem) :
	coupling::solvers::AbstractCouetteSolver<3>(),
	_rank(rank),
	_channelHeight(channelHeight),
	_dx(dx), 
	_dt(dt), 
	_kinVisc(kinVisc), 
	_plotEveryTimestep(plotEveryTimestep), 
	_filestem(filestem),
	_domainSize(std::floor(_channelHeight/dx)),
	_numberOfVertices(std::pow(_domainSize, 3))
	{
		std::cout << "PreciceSolver constructor call" << std::endl;	
		if(skipRank()){return;}
		_interface = new precice::SolverInterface("MAMICO","precice-config.xml",0,1);
		int dim = _interface->getDimensions();
		if (dim!=3)
		{
			std::cout << "main: preCICE dimension should be 3" << std::endl;
			exit(EXIT_FAILURE);
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
		for (int x = 0; x < _domainSize; x++)
		{
			for (int y = 0; y < _domainSize; y++)
			{
				for (int z = 0; z < _domainSize; z++)
				{
					coords[dim*get(x,y,z)] = x * _dx;
					coords[dim*get(x,y,z)+1] = y * _dx;
					coords[dim*get(x,y,z)+2] = z * _dx;			
				}
			}
		}
		_vertexIDs = new int[_numberOfVertices];
		_interface->setMeshVertices(meshID, _numberOfVertices, coords, _vertexIDs); 
		delete[] coords;
		_precice_dt = _interface->initialize();
		std::cout << "PreciceSolver constructor called" << std::endl;
	}

    /** @brief a simple destructor 
    */
    virtual ~PreciceSolver()
    {
		if (_vel !=NULL){delete [] _vel; _vel=NULL;}
		if (_density!=NULL){delete [] _density; _density=NULL;}
		if (_interface!=NULL){_interface->finalize(); delete _interface; _interface=NULL;}
    }
    
    /** @brief changes the velocity at the moving wall (z=0)
     *  @param wallVelocity the velocity will be set at the moving wall 
     */
	void setWallVelocity(const tarch::la::Vector<3,double> wallVelocity) override
	{
		_wallVelocity = wallVelocity;
	}
	
	/** @brief gets the velocity at a given position
	 *  @param pos a position within the continuum domain
	 *  @returns the velocity vector  
	 */
	tarch::la::Vector<3,double> getVelocity(tarch::la::Vector<3,double> pos) const override
	{
		// compute index for respective cell (_dx+... for ghost cells); use coords to store local cell coordinates
		tarch::la::Vector<3,unsigned int> coords;
		for (unsigned int d = 0; d < 3; d++) 
			coords[d] =  (unsigned int) pos[d]/_dx;
		const int index = get(coords[0],coords[1],coords[2]);
		tarch::la::Vector<3,double> vel(0.0);
		// extract velocity
		for (int d = 0; d < 3; d++)
			vel[d] = _vel[3*index+d];
		return vel;
	}
	
	/** @brief returns density at a certain position
	 *  @param pos position for which the density will be returned
	 *  @returns a density vector 
	 */
	double getDensity(tarch::la::Vector<3,double> pos) const
	{
		tarch::la::Vector<3,unsigned int> coords;
		for (unsigned int d = 0; d < 3; d++)
			coords[d] =  (unsigned int) pos[d]/_dx;
		const int index = get(coords[0],coords[1],coords[2]);
		return _density[index];
	}
	
    /** @brief advances one time step dt in time and triggers vtk plot if required 
    */
	void advance(double dt) override 
	{
		if(skipRank()){return;}
		int dim = _interface->getDimensions();
		int meshID = _interface->getMeshID("MamicoMesh");
		int exvsID = _interface->getDataID("ExternalVelocities", meshID); 
		double* velocities = new double[_numberOfVertices*dim];
		_interface->readBlockVectorData(exvsID, _numberOfVertices, _vertexIDs, velocities);
		for (int vertexIndex = 0; vertexIndex < _numberOfVertices; vertexIndex++)
		{
			std::cout << "PreciceSolver::advance: read velocities[" << vertexIndex << "] = [" << velocities[dim*vertexIndex] << "," << velocities[dim*vertexIndex+1] << "," << velocities[dim*vertexIndex+2] << "]" << std::endl;
		}
		std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
		//while (_interface->isCouplingOngoing()) {
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
			std::cout << "PreciceSolver::advance preCICE solver advanced with dt = " << computed_dt << std::endl;
    	//}
	}

  private:    
    /** @brief returns linearized index
     *  @returns the linearized index */
    int get(int x,int y,int z) const
    {
      return x + _domainSize*(y+_domainSize*z);
    }

    /** @brief returns true, if this rank is not of relevance for the LB simulation
     *  @returns a bool, which indicates if the rank shall not do anything (true) or not (false) 
     */
	bool skipRank() const {
    	return !(_rank==0);
  	}

    /** @brief rank of this process */
    const int _rank;
    /** @brief the height and width of the channel in z and y direction */
    const double _channelHeight; //
    /** @brief mesh size, dx=dy=dz */
    const double _dx;
    /** @brief time step*/
    const double _dt;
    /** @brief kinematic viscosity of the fluid */
    const double _kinVisc;
    /** @brief number of time steps between vtk plots */
    const int _plotEveryTimestep;
    /** @brief file stem for vtk plot */
    const std::string _filestem;
    /** @brief domain size in x-y-z-direction */
    const int _domainSize;
     /** @brief number of vertices */
    const int _numberOfVertices;
    /** @brief time step counter */
    int _counter{0};
    /** @brief velocity field */
    double *_vel{NULL};
    /** @brief density field */
    double *_density{NULL};
    /** @brief  velocity of moving wall of Couette flow */
	tarch::la::Vector<3,double> _wallVelocity;
	/** @brief precice solver interface */
    precice::SolverInterface *_interface;
    /** @brief preCICE time step */
	double _precice_dt;
	/** @brief preCICE vertex ids **/
	int* _vertexIDs;
};

#endif// _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
