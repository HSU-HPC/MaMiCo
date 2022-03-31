// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_

#include "tarch/la/Vector.h"
#include <cmath>

namespace coupling {
namespace solvers {
/** interface for Couette solvers. advance(dt) advances the solver in time.
 * getVelocity(pos) returns the velocity at position pos. */
template <unsigned int dim> class AbstractCouetteSolver {
public:
  virtual ~AbstractCouetteSolver() {}
  virtual void advance(double dt) = 0;
  virtual tarch::la::Vector<dim, double> getVelocity(tarch::la::Vector<dim, double> pos) const = 0;
  virtual void setWallVelocity(const tarch::la::Vector<dim, double> wallVelocity) = 0;
};

template <unsigned int dim> class CouetteSolver;
} // namespace solvers
} // namespace coupling

/** implements an analytic Couette flow solver.
 *  In our scenario, the lower wall is accelerated and the upper wall stands
 * still.
 *  The lower wall is located at zero height.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::solvers::CouetteSolver : public coupling::solvers::AbstractCouetteSolver<dim> {
public:
  CouetteSolver(const double &channelheight, const double &wallVelocity, const double kinVisc)
      : AbstractCouetteSolver<dim>(), _channelheight(channelheight), _wallVelocity(wallVelocity), _kinVisc(kinVisc), _time(0.0) {}
  virtual ~CouetteSolver() {}

  /** advances one time step dt in time */
  virtual void advance(double dt) { _time += dt; }

  /** returns the velocity vector at a certain channel height pos[dim-1], given
   * by:
   *  u(z,t)= V(1-z/H) - 2V/pi*sum_k=1^infty 1/k*sin(k*pi*z/H)*exp(-k^2 pi^2/H^2
   * * nu*t)
   */
  virtual tarch::la::Vector<dim, double> getVelocity(tarch::la::Vector<dim, double> pos) const {
    tarch::la::Vector<dim, double> v(0.0);
    v[0] = _wallVelocity * (1.0 - pos[dim - 1] / _channelheight);
    const double pi = 3.141592653589793238;
    double sum = 0.0;
    for (int k = 1; k < 30; k++) {
      sum += 1.0 / k * sin(k * pi * pos[dim - 1] / _channelheight) * exp(-k * k * pi * pi / (_channelheight * _channelheight) * _kinVisc * _time);
    }
    v[0] = v[0] - 2.0 * _wallVelocity / pi * sum;
    return v;
  }

  virtual void setWallVelocity(const tarch::la::Vector<dim, double> wallVelocity) { _wallVelocity = wallVelocity[0]; }

private:
  const double _channelheight; // height of couette channel
  double _wallVelocity;        // velocity of moving wall
  const double _kinVisc;       // kinematic viscosity
  double _time;                // current time
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVER_H_
