// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _DUMMY_SOLVER_H
#define _DUMMY_SOLVER_H

#include "tarch/la/Vector.h"
#include <iostream>


/** dummy solver, implementing a flow solver with parabolic velocity profile and constant density everywhere.
 *  @author Rahul Arora
 */
class DummySolver {
public:
   DummySolver(int nx, int ny, int nz, double vmax):_Nx(nx), _Ny(ny), _Nz(nz), _Vmax(vmax){
	_domain = new double[_Nx*_Ny*_Nz*4];
	for(unsigned int x=0; x<_Nx; x++){
		for(unsigned int y=0; y<_Ny; y++){
			for(unsigned int z=0; z<_Nz; z++){
				_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*0]=0.6;				// density=1.0
				_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*1]=(4.0*_Vmax*((double)(y)/(_Ny-1))*(1-((double)(y)/(_Ny-1))));	// Vx=Vmax*y*(Ny-y)
                                _domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*2]=0.0;		   		// Vy=0.0
				_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*3]=0.0;   			// Vz=0.0
			}
		}
	}
   };

   ~DummySolver(){
	delete[] _domain;
   }

   void setVelocity(tarch::la::Vector<3,double> velocity, int x, int y, int z){
        _domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*1]=velocity[0];       
        _domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*2]=velocity[1];
        _domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*3]=velocity[2];
   }

   tarch::la::Vector<3,double> getVelocity(int x, int y, int z) const{
   	tarch::la::Vector<3,double> velocity(0.0);
	velocity[0]=_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*1];
        velocity[1]=_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*2];
        velocity[2]=_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*3];
	return velocity;
   }

   void setDensity(double density, int x, int y, int z){
	_domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*0]=density;
   }

   double getDensity(int x, int y, int z){
	return _domain[x+_Nx*y+_Nx*_Ny*z+_Nx*_Ny*_Nz*0];
   }

private: 
   unsigned int _Nx, _Ny, _Nz;
   double _Vmax;
   double* _domain;
};

#endif //_DUMMY_SOLVER_H
