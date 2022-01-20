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
	  const double channelheight,
	  const double dx,
	  const double dt,
	  const double kinVisc,
	  const int plotEveryTimestep,
	  const std::string filestem) :
	coupling::solvers::AbstractCouetteSolver<3>(),
	_rank(rank),
	_channelheight(channelheight),
	_dx(dx), 
	_dt(dt), 
	_kinVisc(kinVisc), 
	_plotEveryTimestep(plotEveryTimestep), 
	_filestem(filestem)
	{
		std::cout << "PreciceSolver constructor call" << std::endl;	
		if(skipRank()){return;}
		
		_vel = new double[3*(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
		_density = new double[(_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2)];
		// zero velocity, unit density;
		for (int i = 0; i < (_domainSizeX+2)*(_domainSizeY+2)*(_domainSizeZ+2); i++)
		{
			for (int d = 0; d <  3; d++)
				_vel[i*3+d]  = (double)0.0;
			_density[i] = 1.0;
		}
		_interface = new precice::SolverInterface("MAMICO","precice-config.xml",0,1);
		int dim = _interface->getDimensions();
		int meshID = _interface->getMeshID("MamicoMesh");
		int vertexSize = 10;
		// determine vertexSize
		double* coords = new double[vertexSize*dim]; // coords of coupling vertices 
		// determine coordinates
		int* vertexIDs = new int[vertexSize];
		_interface->setMeshVertices(meshID, vertexSize, coords, vertexIDs); 
		delete[] coords;
		delete[] vertexIDs;

		/*int mvsID = precice->getDataID("MamicoVelocities", meshID); 
		int dvsID = precice->getDataID("DummmyVelocities", meshID); 
		double* mvs = new double[vertexSize*dim];
		double* dvs = new double[vertexSize*dim];*/

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
			coords[d] =  (unsigned int) ((_dx+pos[d])/_dx);
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
			coords[d] =  (unsigned int) ((_dx+pos[d])/_dx);
		const int index = get(coords[0],coords[1],coords[2]);
		return _density[index];
	}
	
    /** @brief advances one time step dt in time and triggers vtk plot if required 
    */
	void advance(double dt) override 
	{
		/*
		if (skipRank()){return;}
		const int timesteps=floor( dt/_dt+0.5 );
		if ( fabs(timesteps*_dt-dt)/_dt > 1.0e-8 )
		{
			std::cout << "ERROR PreCICESolver::advance(): time steps and dt do not match!" << std::endl; 
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i < timesteps; i++)
		{
			setBeyondWall();
			update();
			plot();
			_counter++;
		}
		*/
		if(skipRank()){return;}
		std::cout << "PreciceSolver::advance called with dt = " << dt << std::endl;
		//while (_interface->isCouplingOngoing()) {
			double computed_dt = std::min(_precice_dt, dt);
			_precice_dt = _interface->advance(computed_dt);
			std::cout << "PreciceSolver::advance preCICE solver advanced with dt = " << computed_dt << std::endl;
    	//}
	}

  private:
    /** @brief determines the local domain size on this rank where channelheight is the domain length in direction d.
     *  @returns the size of the domain */
    int getDomainSize(double channelheight, double dx) const 
    {
      return floor( (channelheight+0.5)/dx );
    }
    
    /** @brief returns linearized index and performs checks in debug mode
     *  @returns the linearized index */
    int get(int x,int y,int z) const
    {
      return x + (_domainSizeX+2)*(y+(_domainSizeY+2)*z);
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
    const double _channelheight; //
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
    /** @brief domain size in x-direction */
    const int _domainSizeX { getDomainSize(_channelheight,_dx) };
    /** @brief domain size in y-direction */
    const int _domainSizeY { getDomainSize(_channelheight,_dx) };
    /** @brief domain size in z-direction */
    const int _domainSizeZ { getDomainSize(_channelheight,_dx) };
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
};

#endif// _MOLECULARDYNAMICS_COUPLING_SOLVERS_PRECICESOLVER_H_
