#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_FINITEDIFFERENCE_H_

#include "coupling/solvers/NumericalSolver.h"
#include "tarch/la/Vector.h"
#include <cmath>

namespace coupling {
namespace solvers {
class FiniteDifferenceSolver;
}
} // namespace coupling

/** In our scenario, the lower wall is accelerated and the upper wall stands
 * still. The lower wall is located at zero height. The grid is just a simple
 * cubic, equidistant mesh.
 *  @brief implements a simple one-dimensional finite-diffference solver for the
 * Couette flow.
 *  @author Helene Wittenberg  */
class coupling::solvers::FiniteDifferenceSolver : public coupling::solvers::NumericalSolver {
public:
  /** @brief a simple constructor
   *  @param channelheight the width and height of the channel in y and z
   * direction
   *  @param kinVisc the kinematic viscosity of the fluid
   *  @param wallVelocity velocity at the moving wall, refers to Couette
   * scenario
   *  @param dx the spacial step size, and equidistant grid is applied
   *  @param dt the time step
   *  @param plotEveryTimestep the time step interval for plotting data;
   *                           4 means, every 4th time step is plotted
   *  @param filestem the name of the plotted file
   *  @param processes defines on how many processes the solver will run;
   *                   1,1,1 - sequential run - 1,2,2 = 1*2*2 = 4 processes */
  FiniteDifferenceSolver(const double channelheight, tarch::la::Vector<3, double> wallVelocity, const double kinVisc, const double dx, const double dt,
                         const int plotEveryTimestep, const std::string filestem, const tarch::la::Vector<3, unsigned int> processes)
      : coupling::solvers::NumericalSolver(channelheight, dx, dt, kinVisc, plotEveryTimestep, filestem, processes), _omega(dt * kinVisc / (dx * dx)),
        _wallVelocity(wallVelocity) {
    if (_processes[0] > 1 || _processes[1] > 1 || _processes[2] > 1) {
      _processes[0] = 1;
      _processes[1] = 1;
      _processes[2] = 1;
      std::cout << "The FiniteDifferenceSolver was requested in a parallel manner. It will be run sequentially, since it doesn't support parallel runs. "
                << std::endl;
    }
    // return if required
    if (skipRank()) {
      return;
    }
    _velold = new double[3 * (_domainSizeX + 2) * (_domainSizeY + 2) * (_domainSizeZ + 2)];
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Domain size=" << _domainSizeX << "," << _domainSizeY << "," << _domainSizeZ << std::endl;
    std::cout << "tau=" << 1.0 / _omega << std::endl;
    std::cout << "wallVelocity=" << _wallVelocity << std::endl;
    for (int z = 0; z < _domainSizeZ + 2; z++) {
      for (int y = 0; y < _domainSizeY + 2; y++) {
        for (int x = 0; x < _domainSizeX + 2; x++) {
          std::cout << x << "," << y << "," << z << "FLAG=" << _flag[get(x, y, z)] << std::endl;
        }
      }
    }
#endif
    // check pointers
    if ((!_vel) || (!_density) || (!_flag)) {
      std::cout << "ERROR FiniteDifferenceSolver: nullptr!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** @brief a simple destructor  */
  ~FiniteDifferenceSolver() {
    if (_velold) {
      delete[] _velold;
      _velold = nullptr;
    }
  }

  /** @brief advances one time step dt in time and triggers vtk plot if required
   */
  void advance(double dt) override {
    if (skipRank()) {
      return;
    }
    const int timesteps = floor(dt / _dt + 0.5);
    if (fabs(timesteps * _dt - dt) / _dt > 1.0e-8) {
      std::cout << "ERROR FiniteDifferenceSolver::advance(): time steps and dt do not match!" << std::endl;
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < timesteps; i++) {
      setBeyondWall();
      update();
      plot();
      _counter++;
    }
  }

  /** @brief gets the velocity at a given position
   *  @param pos a position within the continuum domain
   *  @returns the velocity vector  */
  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const override {
    const tarch::la::Vector<3, double> domainOffset(_coords[0] * _dx * _avgDomainSizeX, _coords[1] * _dx * _avgDomainSizeY, _coords[2] * _dx * _avgDomainSizeZ);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if ((pos[0] < domainOffset[0]) || (pos[0] > domainOffset[0] + _domainSizeX * _dx) || (pos[1] < domainOffset[1]) ||
        (pos[1] > domainOffset[1] + _domainSizeY * _dx) || (pos[2] < domainOffset[2]) || (pos[2] > domainOffset[2] + _domainSizeZ * _dx)) {
      std::cout << "ERROR FiniteDifferenceSolver::getVelocity(): Position " << pos << " out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    // compute index for respective cell (_dx+... for ghost cells); use coords
    // to store local cell coordinates
    tarch::la::Vector<3, unsigned int> coords;
    for (unsigned int d = 0; d < 3; d++) {
      coords[d] = (unsigned int)((_dx + pos[d] - domainOffset[d]) / _dx);
    }
    const int index = get(coords[0], coords[1], coords[2]);
    tarch::la::Vector<3, double> vel(0.0);
    // extract velocity
    for (int d = 0; d < 3; d++) {
      vel[d] = _vel[3 * index + d];
    }
    return vel;
  }

  /** @brief changes the velocity at the moving wall (z=0)
   *  @param wallVelocity the velocity will be set at the moving wall */
  virtual void setWallVelocity(const tarch::la::Vector<3, double> wallVelocity) override { _wallVelocity = wallVelocity; }

  /** @brief returns density at a certain position
   *  @param pos position for which the density will be returned
   *  @returns a density vector */
  double getDensity(tarch::la::Vector<3, double> pos) const override {
    tarch::la::Vector<3, unsigned int> coords;
    const tarch::la::Vector<3, double> domainOffset(_coords[0] * _dx * _avgDomainSizeX, _coords[1] * _dx * _avgDomainSizeY, _coords[2] * _dx * _avgDomainSizeZ);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if ((pos[0] < domainOffset[0]) || (pos[0] > domainOffset[0] + _domainSizeX * _dx) || (pos[1] < domainOffset[1]) ||
        (pos[1] > domainOffset[1] + _domainSizeY * _dx) || (pos[2] < domainOffset[2]) || (pos[2] > domainOffset[2] + _domainSizeZ * _dx)) {
      std::cout << "ERROR FiniteDifferenceSolver::getDensity(): Position " << pos << " out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    for (unsigned int d = 0; d < 3; d++) {
      coords[d] = (unsigned int)((_dx + pos[d] - domainOffset[d]) / _dx);
    }
    const int index = get(coords[0], coords[1], coords[2]);
    return _density[index];
  }

  /** @brief applies the values received from the MD-solver within the
   * conntinuum solver
   *  @param recvBuffer holds the data from the md solver
   *  @param recvIndice the indices to connect the data from the buffer with
   * coupling cells */

  void setMDBoundaryValues(coupling::datastructures::FlexibleCellContainer<3>& md2macroBuffer) override {
    if (skipRank()) {
      return;
    }
    for (auto pair : md2macroBuffer) {
      I01 idx;
      coupling::datastructures::CouplingCell<3>* couplingCell;
      std::tie(couplingCell, idx) = pair;
      // determine cell index of this cell in continuum domain
      tarch::la::Vector<3, unsigned int> globalCellCoords{I01{idx}.get()};
      globalCellCoords[0] = (globalCellCoords[0] + _offset[0]) - _coords[0] * _avgDomainSizeX;
      globalCellCoords[1] = (globalCellCoords[1] + _offset[1]) - _coords[1] * _avgDomainSizeY;
      globalCellCoords[2] = (globalCellCoords[2] + _offset[2]) - _coords[2] * _avgDomainSizeZ;
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "Process coords: " << _coords << ":  GlobalCellCoords for index ";
      std::cout << I01{idx} << ": " << globalCellCoords << std::endl;
#endif
      const int index = get(globalCellCoords[0], globalCellCoords[1], globalCellCoords[2]);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (_flag[index] != MD_BOUNDARY) {
        std::cout << "ERROR FiniteDifferenceSolver::setMDBoundaryValues(): Cell " << index << " is no MD boundary cell!" << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      // set velocity value in MD boundary cell (before streaming); the boundary velocities are interpolated between the neighbouring and this cell. This
      // interpolation is valid for FLUID-MD_BOUNDARY neighbouring relations only. determine local velocity received from MaMiCo and convert it to LB units;
      // store the velocity in _vel
      tarch::la::Vector<3, double> localVel((1.0 / couplingCell->getMacroscopicMass()) * couplingCell->getMacroscopicMomentum());
      for (unsigned int d = 0; d < 3; d++) {
        _vel[3 * index + d] = localVel[d];
      }
    }
  }

private:
  /** there is no grid point on the wall, since the wall is located directly
   * between the cells on the boundary. Therefore so set the boundary condition,
   * the velocity boundary condition is applied by setting the velocity in the
   * outer cell
   *  @brief updates the velocities at boundaries*/
  void setBeyondWall() {
#pragma omp parallel for
    for (int i = 0; i < _yO; i++) {
      _vel[3 * i] = 2 * _wallVelocity[0] - _vel[3 * (i + _yO)];
    }
  }

  /** @brief update the velocity in the cells which are flagged as FLUID */
  void update() {
    double* swap = _velold;
    _velold = _vel;
    _vel = swap;
#pragma omp parallel for
    for (int i = _yO; i < (_domainSizeX + 2) * (_domainSizeY + 2) * (_domainSizeZ + 2) - _yO; i++) {
      if (_flag[i] == FLUID) {
        _vel[3 * i] = _omega * (_velold[3 * (i - _yO)] - 2 * _velold[3 * i] + _velold[3 * (i + _yO)]) + _velold[3 * i];
      }
    }
  }

  /** @brief factor for the finite difference stencil = dt*kinVisc/(dx*dx) */
  const double _omega;
  /** @brief  velocity of moving wall of Couette flow */
  tarch::la::Vector<3, double> _wallVelocity;
  /** @brief the velocity field from the last time step  */
  double* _velold{nullptr};
};
#endif
