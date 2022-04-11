// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_

#include "tarch/la/Vector.h"
#include <cmath>

namespace coupling {
  /** @brief all numerical solvers are defined in the namespace, and their interfaces  */
  namespace solvers{
    /** @brief interface for continuum/macro fluid solvers for the Couette scenario
     *  @author Philipp Neumann
     *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2, or 3*/
    template<unsigned int dim>
    class AbstractCouetteSolver {
      public:
        /** @brief a dummy destructor
         *  @todo Why a dummy destructor buut not a constructor */
        virtual ~AbstractCouetteSolver(){}
        /** @brief advances the solver in time
         *  @param dt the solver will be advanced from the current time t to t+dt */
        virtual void advance(double dt)=0;
        /** @brief returns the current velocity at the given position
         *  @param pos position to return the velocity for
         *  @returns a velocity vector */
        virtual tarch::la::Vector<dim,double> getVelocity(tarch::la::Vector<dim,double> pos) const = 0;
        /** @brief changes the velocity at the moving for, refers to Couette scenario
         *  @param wallVelocity value to which the veloctiy will be set */
        virtual void setWallVelocity(const tarch::la::Vector<dim,double> wallVelocity);
        virtual void writeError(tarch::la::Vector<3,double>& lowerLeftCorner, tarch::la::Vector<3,double>& upperRightCorner)const{}
    };

    template<unsigned int dim>
    class CouetteSolver;
  }
}

/** In our scenario, the lower wall is accelerated and the upper wall stands still.
 *  The lower wall is located at zero height.
 *  @brief implements an analytic Couette flow solver.
 *  @author Philipp Neumann
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2, or 3 */
template<unsigned int dim>
class coupling::solvers::CouetteSolver: public coupling::solvers::AbstractCouetteSolver<dim> {
  public:
    /** @brief a simple constructor
     *  @param channelheight the height and width of the channel in y and z direction
     *  @param wallVelocity the velocity at the moving wall
     *  @param kinVisc the kinematic viscosity of the fluid  */
    CouetteSolver(
      const double &channelheight,
      const double &wallVelocity,
      const double kinVisc
    ): AbstractCouetteSolver<dim>(), _channelheight(channelheight), _wallVelocity(wallVelocity), _kinVisc(kinVisc), _time(0.0){}

    /** @brief a dummy destructor */
    virtual ~CouetteSolver(){}

    /** @brief advances one time step dt in time
     *  @param dt size of the time step to advance */
    virtual void advance(double dt){ _time += dt; }

    /** for the first entry of the vector is the analytic solution of the Couette scenario returned, given by:
     *  u(z,t)= V(1-z/H) - 2V/pi*sum_k=1^infty 1/k*sin(k*pi*z/H)*exp(-k^2 pi^2/H^2 * nu*t)
     *  The other two entries are 0
     *  @brief returns the velocity vector at a certain channel position
     *  @param pos the position to return the velocity for
     *  @returns a velocity vector */
    virtual tarch::la::Vector<dim,double> getVelocity(tarch::la::Vector<dim,double> pos) const {
      tarch::la::Vector<dim,double> v(0.0);
      v[0] = _wallVelocity*(1.0-pos[dim-1]/_channelheight);
      const double pi = 3.141592653589793238;
      double sum = 0.0;
      for (int k = 1; k < 30; k++){
        sum += 1.0/k * sin(k*pi*pos[dim-1]/_channelheight) * exp(-k*k * pi*pi/(_channelheight*_channelheight) * _kinVisc * _time);
      }
      v[0] = v[0] - 2.0*_wallVelocity/pi * sum;
      return v;
    }

    /** @brief changes the velocity at the moving for, refers to Couette scenario
     *  @param wallVelocity value to which the veloctiy will be set */
    virtual void setWallVelocity(const tarch::la::Vector<dim,double> wallVelocity){
      _wallVelocity = wallVelocity[0];
    }

  private:
    /** @brief height of couette channel */
    const double _channelheight;
    /** @brief velocity of moving wall */
    double _wallVelocity;
    /** @brief kinematic viscosity */
    const double _kinVisc;
    /** @brief current time */
    double _time;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_
