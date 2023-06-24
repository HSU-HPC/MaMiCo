// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "coupling/solvers/NumericalSolver.h"

namespace coupling {
namespace solvers {
class LBCouetteSolver;
}
} // namespace coupling

/** In our scenario, the lower wall is accelerated and the upper wall stands
 * still. The lower wall is located at zero height.
 *  @brief implements a three-dimensional Lattice-Boltzmann Couette flow solver.
 *  @author Philipp Neumann  */
class coupling::solvers::LBCouetteSolver : public coupling::solvers::NumericalSolver {
public:
  /** @brief a simple constructor
   *  @param channelheight the width and height of the channel in y and z
   * direction
   *  @param wallVelocity velocity at the moving wall, refers to Couette
   * scenario
   *  @param dx the spacial step size, and equidistant grid is applied
   *  @param dt the time step
   *  @param kinVisc the kinematic viscosity of the fluid
   *  @param plotEveryTimestep the time step interval for plotting data;
   *                           4 means, every 4th time step is plotted
   *  @param filestem the name of the plotted file
   *  @param processes defines on how many processes the solver will run;
   *                   1,1,1 - sequential run - 1,2,2 = 1*2*2 = 4 processes
   *  @param numThreads number of OpenMP threads */
  LBCouetteSolver(const double channelheight, tarch::la::Vector<3, double> wallVelocity, const double kinVisc, const double dx, const double dt,
                  const int plotEveryTimestep, const std::string filestem, const tarch::la::Vector<3, unsigned int> processes,
                  const unsigned int numThreads = 1)
      : coupling::solvers::NumericalSolver(channelheight, dx, dt, kinVisc, plotEveryTimestep, filestem, processes),
        _omega(1.0 / (3.0 * (kinVisc * dt / (dx * dx)) + 0.5)), _wallVelocity((dt / dx) * wallVelocity) {
    // return if required
    if (skipRank()) {
      return;
    }
    _pdf1 = new double[19 * (_domainSizeX + 2) * (_domainSizeY + 2) * (_domainSizeZ + 2)];
    _pdf2 = new double[19 * (_domainSizeX + 2) * (_domainSizeY + 2) * (_domainSizeZ + 2)];
#if defined(_OPENMP)
    omp_set_num_threads(numThreads);
#endif
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
    if ((_pdf1 == NULL) || (_pdf2 == NULL) || (_vel == NULL) || (_density == NULL) || (_flag == NULL)) {
      std::cout << "ERROR LBCouetteSolver: NULL ptr!" << std::endl;
      exit(EXIT_FAILURE);
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    if ((_sendBufferX == NULL) || (_recvBufferX == NULL) || (_sendBufferY == NULL) || (_recvBufferY == NULL) || (_sendBufferZ == NULL) ||
        (_recvBufferZ == NULL)) {
      std::cout << "ERROR LBCouetteSolver: NULL ptr in send/recv!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
// init everything with lattice weights
#pragma omp parallel for
    for (int i = 0; i < (_domainSizeX + 2) * (_domainSizeY + 2) * (_domainSizeZ + 2); i++) {
      for (int q = 0; q < 19; q++) {
        _pdf1[get(i) * 19 + q] = _W[q];
        _pdf2[get(i) * 19 + q] = _W[q];
      }
    }
  }

  /** @brief a simple destructor */
  virtual ~LBCouetteSolver() {
    if (_pdf1 != NULL) {
      delete[] _pdf1;
      _pdf1 = NULL;
    }
    if (_pdf2 != NULL) {
      delete[] _pdf2;
      _pdf2 = NULL;
    }
    if (_vel != NULL) {
      delete[] _vel;
      _vel = NULL;
    }
    if (_density != NULL) {
      delete[] _density;
      _density = NULL;
    }
    if (_flag != NULL) {
      delete[] _flag;
      _flag = NULL;
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    if (_sendBufferX != NULL) {
      delete[] _sendBufferX;
      _sendBufferX = NULL;
    }
    if (_sendBufferY != NULL) {
      delete[] _sendBufferY;
      _sendBufferY = NULL;
    }
    if (_sendBufferZ != NULL) {
      delete[] _sendBufferZ;
      _sendBufferZ = NULL;
    }
    if (_recvBufferX != NULL) {
      delete[] _recvBufferX;
      _recvBufferX = NULL;
    }
    if (_recvBufferY != NULL) {
      delete[] _recvBufferY;
      _recvBufferY = NULL;
    }
    if (_recvBufferZ != NULL) {
      delete[] _recvBufferZ;
      _recvBufferZ = NULL;
    }
#endif
  }

  /** @brief advances one time step dt in time and triggers vtk plot if required
   */
  void advance(double dt) override {
    if (skipRank()) {
      return;
    }
    const int timesteps = floor(dt / _dt + 0.5);
    if (fabs(timesteps * _dt - dt) / _dt > 1.0e-8) {
      std::cout << "ERROR LBCouetteSolver::advance(): time steps and dt do not match!" << std::endl;
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < timesteps; i++) {
      plot();
      collidestream();
      communicate(); // exchange between neighbouring MPI subdomains
      _counter++;
    }
  }

  /** @brief applies the values received from the MD-solver within the
   * conntinuum solver
   *  @param recvBuffer holds the data from the md solver
   *  @param recvIndice the indices to connect the data from the buffer with
   * macroscopic cells
   *  @param indexConversion instance of the indexConversion */
  void setMDBoundaryValues(std::vector<coupling::datastructures::MacroscopicCell<3>*>& recvBuffer, const unsigned int* const recvIndices,
                           const coupling::IndexConversion<3>& indexConversion) override {
    if (skipRank()) {
      return;
    }
    // loop over all received cells
    const unsigned int size = (unsigned int)recvBuffer.size();
    for (unsigned int i = 0; i < size; i++) {
      // determine cell index of this cell in LB domain
      tarch::la::Vector<3, unsigned int> globalCellCoords = indexConversion.getGlobalVectorCellIndex(recvIndices[i]);
      globalCellCoords[0] = (globalCellCoords[0] + _offset[0]) - _coords[0] * _avgDomainSizeX;
      globalCellCoords[1] = (globalCellCoords[1] + _offset[1]) - _coords[1] * _avgDomainSizeY;
      globalCellCoords[2] = (globalCellCoords[2] + _offset[2]) - _coords[2] * _avgDomainSizeZ;
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "Process coords: " << _coords << ":  GlobalCellCoords for index " << indexConversion.getGlobalVectorCellIndex(recvIndices[i]) << ": "
                << globalCellCoords << std::endl;
#endif
      const int index = get(globalCellCoords[0], globalCellCoords[1], globalCellCoords[2]);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      if (_flag[index] != MD_BOUNDARY) {
        std::cout << "ERROR LBCouetteSolver::setMDBoundaryValues(): Cell " << index << " is no MD boundary cell!" << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      // set velocity value and pdfs in MD boundary cell (before streaming); the
      // boundary velocities are interpolated between the neighbouring and this
      // cell. This interpolation is valid for FLUID-MD_BOUNDARY neighbouring
      // relations only. determine local velocity received from MaMiCo and
      // convert it to LB units; store the velocity in _vel
      tarch::la::Vector<3, double> localVel((1.0 / recvBuffer[i]->getMacroscopicMass()) * (_dt / _dx) * recvBuffer[i]->getMacroscopicMomentum());
      for (unsigned int d = 0; d < 3; d++) {
        _vel[3 * index + d] = localVel[d];
      }
      // loop over all pdfs and set them according to interpolated moving-wall
      // conditions
      for (unsigned int q = 0; q < 19; q++) {
        // index of neighbour cell; only if cell is located inside local domain
        if (((int)globalCellCoords[0] + _C[q][0] > 0) && ((int)globalCellCoords[0] + _C[q][0] < _domainSizeX + 1) &&
            ((int)globalCellCoords[1] + _C[q][1] > 0) && ((int)globalCellCoords[1] + _C[q][1] < _domainSizeY + 1) &&
            ((int)globalCellCoords[2] + _C[q][2] > 0) && ((int)globalCellCoords[2] + _C[q][2] < _domainSizeZ + 1)) {
          const int nbIndex = get((_C[q][0] + globalCellCoords[0]), (_C[q][1] + globalCellCoords[1]), (_C[q][2] + globalCellCoords[2]));
          const tarch::la::Vector<3, double> interpolVel(0.5 * (_vel[3 * index] + _vel[3 * nbIndex]), 0.5 * (_vel[3 * index + 1] + _vel[3 * nbIndex + 1]),
                                                         0.5 * (_vel[3 * index + 2] + _vel[3 * nbIndex + 2]));
          _pdf1[19 * index + q] =
              _pdf1[19 * nbIndex + 18 - q] -
              6.0 * _W[q] * _density[nbIndex] * (_C[18 - q][0] * interpolVel[0] + _C[18 - q][1] * interpolVel[1] + _C[18 - q][2] * interpolVel[2]);
        }
      }
    }
  }

  /** @brief returns velocity at a certain position
   *  @param pos position for which the velocity will be returned
   *  @returns the velocity vector for the position */
  tarch::la::Vector<3, double> getVelocity(tarch::la::Vector<3, double> pos) const override {
    tarch::la::Vector<3, unsigned int> coords;
    const tarch::la::Vector<3, double> domainOffset(_coords[0] * _dx * _avgDomainSizeX, _coords[1] * _dx * _avgDomainSizeY, _coords[2] * _dx * _avgDomainSizeZ);
    // check pos-data for process locality (todo: put this in debug mode in
    // future releases)
    if ((pos[0] < domainOffset[0]) || (pos[0] > domainOffset[0] + _domainSizeX * _dx) || (pos[1] < domainOffset[1]) ||
        (pos[1] > domainOffset[1] + _domainSizeY * _dx) || (pos[2] < domainOffset[2]) || (pos[2] > domainOffset[2] + _domainSizeZ * _dx)) {
      std::cout << "ERROR LBCouetteSolver::getVelocity(): Position " << pos << " out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // compute index for respective cell (_dx+... for ghost cells); use coords
    // to store local cell coordinates
    for (unsigned int d = 0; d < 3; d++) {
      coords[d] = (unsigned int)((_dx + pos[d] - domainOffset[d]) / _dx);
    }
    const int index = get(coords[0], coords[1], coords[2]);
    tarch::la::Vector<3, double> vel(0.0);
    // extract and scale velocity to "real"=MD units
    for (int d = 0; d < 3; d++) {
      vel[d] = _dx / _dt * _vel[3 * index + d];
    }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Position " << pos << " corresponds to cell: " << coords << "; vel=" << vel << std::endl;
#endif
    return vel;
  }

  /** @brief returns density at a certain position
   *  @param pos position for which the density will be returned
   *  @returns the density vector for the position */
  double getDensity(tarch::la::Vector<3, double> pos) const override {
    tarch::la::Vector<3, unsigned int> coords;
    const tarch::la::Vector<3, double> domainOffset(_coords[0] * _dx * _avgDomainSizeX, _coords[1] * _dx * _avgDomainSizeY, _coords[2] * _dx * _avgDomainSizeZ);
    // check pos-data for process locality (todo: put this in debug mode in
    // future releases)
    if ((pos[0] < domainOffset[0]) || (pos[0] > domainOffset[0] + _domainSizeX * _dx) || (pos[1] < domainOffset[1]) ||
        (pos[1] > domainOffset[1] + _domainSizeY * _dx) || (pos[2] < domainOffset[2]) || (pos[2] > domainOffset[2] + _domainSizeZ * _dx)) {
      std::cout << "ERROR LBCouetteSolver::getDensity(): Position " << pos << " out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // compute index for respective cell (_dx+... for ghost cells); use coords
    // to store local cell coordinates
    for (unsigned int d = 0; d < 3; d++) {
      coords[d] = (unsigned int)((_dx + pos[d] - domainOffset[d]) / _dx);
    }
    const int index = get(coords[0], coords[1], coords[2]);
    return _density[index];
  }

  /** @brief changes the velocity at the moving wall (z=0)
   *  @param wallVelocity the velocity will be set at the moving wall */
  virtual void setWallVelocity(const tarch::la::Vector<3, double> wallVelocity) override { _wallVelocity = (_dt / _dx) * wallVelocity; }

private:
  /** calls stream() and collide() and swaps the fields
   *  @brief collide-stream algorithm for the Lattice-Boltzmann method  */
  void collidestream() {
#pragma omp parallel for
    for (int z = 1; z < _domainSizeZ + 1; z++) {
      for (int y = 1; y < _domainSizeY + 1; y++) {
        for (int x = 1; x < _domainSizeX + 1; x++) {
          const int index = get(x, y, z);
          if (_flag[index] == FLUID) {
            stream(index);
            collide(index, x, y, z);
          }
        }
      }
    }
    // swap fields
    double* swap = _pdf1;
    _pdf1 = _pdf2;
    _pdf2 = swap;
  }

  /** @brief the stream part of the LB algorithm (from pdf1 to pdf2) */
  void stream(int index) {
    const int pI = 19 * index;
    for (int q = 0; q < 9; q++) {
      const int nb = 19 * (_C[q][0] + _C[q][1] * _xO + _C[q][2] * _yO);
      _pdf2[pI + q] = _pdf1[pI + q - nb];
      _pdf2[pI + 18 - q] = _pdf1[pI + 18 - q + nb];
    }
    _pdf2[pI + 9] = _pdf1[pI + 9];
  }

  /** @brieff the collide step within pdf2 */
  void collide(int index, int x, int y, int z) {
    // index of start of cell-local pdfs in AoS
    const int pI = 19 * index;
    // compute and store density, velocity
    double* vel = &_vel[3 * index];
    computeDensityAndVelocity(vel, _density[index], &_pdf2[pI]);
    // collide (BGK); always handle pdfs no. q and inv(q)=18-q in one step
    const double u2 = 1.0 - 1.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    // pdf 0,18
    double cu = -vel[1] - vel[2];
    int nb = -_xO - _yO;
    double feq = _W[0] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI] -= _omega * (_pdf2[pI] - feq);
    boundary(_pdf2, pI, x, y, z, 0, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[0] * _density[index] * cu;
    _pdf2[pI + 18] -= _omega * (_pdf2[pI + 18] - feq);
    boundary(_pdf2, pI, x, y, z, 18, _flag[index - nb], pI - 19 * nb);
    // pdf 1,17
    cu = -vel[0] - vel[2];
    nb = -1 - _yO;
    feq = _W[1] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 1] -= _omega * (_pdf2[pI + 1] - feq);
    boundary(_pdf2, pI, x, y, z, 1, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[1] * _density[index] * cu;
    _pdf2[pI + 17] -= _omega * (_pdf2[pI + 17] - feq);
    boundary(_pdf2, pI, x, y, z, 17, _flag[index - nb], pI - 19 * nb);
    // pdf 2,16
    cu = -vel[2];
    nb = -_yO;
    feq = _W[2] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 2] -= _omega * (_pdf2[pI + 2] - feq);
    boundary(_pdf2, pI, x, y, z, 2, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[2] * _density[index] * cu;
    _pdf2[pI + 16] -= _omega * (_pdf2[pI + 16] - feq);
    boundary(_pdf2, pI, x, y, z, 16, _flag[index - nb], pI - 19 * nb);
    // pdf 3,15
    cu = vel[0] - vel[2];
    nb = 1 - _yO;
    feq = _W[3] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 3] -= _omega * (_pdf2[pI + 3] - feq);
    boundary(_pdf2, pI, x, y, z, 3, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[3] * _density[index] * cu;
    _pdf2[pI + 15] -= _omega * (_pdf2[pI + 15] - feq);
    boundary(_pdf2, pI, x, y, z, 15, _flag[index - nb], pI - 19 * nb);
    // pdf 4,14
    cu = vel[1] - vel[2];
    nb = _xO - _yO;
    feq = _W[4] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 4] -= _omega * (_pdf2[pI + 4] - feq);
    boundary(_pdf2, pI, x, y, z, 4, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[4] * _density[index] * cu;
    _pdf2[pI + 14] -= _omega * (_pdf2[pI + 14] - feq);
    boundary(_pdf2, pI, x, y, z, 14, _flag[index - nb], pI - 19 * nb);
    // pdf 5,13
    cu = -vel[0] - vel[1];
    nb = -1 - _xO;
    feq = _W[5] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 5] -= _omega * (_pdf2[pI + 5] - feq);
    boundary(_pdf2, pI, x, y, z, 5, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[5] * _density[index] * cu;
    _pdf2[pI + 13] -= _omega * (_pdf2[pI + 13] - feq);
    boundary(_pdf2, pI, x, y, z, 13, _flag[index - nb], pI - 19 * nb);
    // pdf 6,12
    cu = -vel[1];
    nb = -_xO;
    feq = _W[6] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 6] -= _omega * (_pdf2[pI + 6] - feq);
    boundary(_pdf2, pI, x, y, z, 6, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[6] * _density[index] * cu;
    _pdf2[pI + 12] -= _omega * (_pdf2[pI + 12] - feq);
    boundary(_pdf2, pI, x, y, z, 12, _flag[index - nb], pI - 19 * nb);
    // pdf 7,11
    cu = vel[0] - vel[1];
    nb = 1 - _xO;
    feq = _W[7] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 7] -= _omega * (_pdf2[pI + 7] - feq);
    boundary(_pdf2, pI, x, y, z, 7, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[7] * _density[index] * cu;
    _pdf2[pI + 11] -= _omega * (_pdf2[pI + 11] - feq);
    boundary(_pdf2, pI, x, y, z, 11, _flag[index - nb], pI - 19 * nb);
    // pdf 8,10
    cu = -vel[0];
    nb = -1;
    feq = _W[8] * _density[index] * (u2 + 3.0 * cu + 4.5 * cu * cu);
    _pdf2[pI + 8] -= _omega * (_pdf2[pI + 8] - feq);
    boundary(_pdf2, pI, x, y, z, 8, _flag[index + nb], pI + 19 * nb);
    feq -= 6.0 * _W[8] * _density[index] * cu;
    _pdf2[pI + 10] -= _omega * (_pdf2[pI + 10] - feq);
    boundary(_pdf2, pI, x, y, z, 10, _flag[index - nb], pI - 19 * nb);
    // pdf 9
    _pdf2[pI + 9] -= _omega * (_pdf2[pI + 9] - _W[9] * _density[index] * u2);
  }

  /** @brief takes care of the correct boundary treatment for the LB method
   *  @param pdf particle distribution function
   *  @param index start index for current cell in pdf-array
   *  @param x the position in x direction of the cell
   *  @param y the position in y direction of the cell
   *  @param z the position in z direction of the cell
   *  @param q distribution function number
   *  @param flag boundary flag of neighbouring cell
   *  @param nbIndex index of neighbouring cell */
  void boundary(double* const pdf, int index, int x, int y, int z, int q, const Flag& flag, int nbIndex) {
    if (flag != FLUID) {
      if (flag == NO_SLIP) {
        // half-way bounce back
        pdf[nbIndex + 18 - q] = pdf[index + q];
      } else if (flag == MOVING_WALL) {
        // half-way bounce back + moving wall acceleration (only x-direction for
        // wall supported at the moment)
        pdf[nbIndex + 18 - q] =
            pdf[index + q] - 6.0 * _W[q] * _density[index / 19] * (_C[q][0] * _wallVelocity[0] + _C[q][1] * _wallVelocity[1] + _C[q][2] * _wallVelocity[2]);
      } else if (flag == PERIODIC) {
        // periodic treatment
        int target[3] = {x, y, z};
        if (target[0] + _C[q][0] == 0) {
          target[0] = _domainSizeX + 1;
        } else if (target[0] + _C[q][0] == _domainSizeX + 1) {
          target[0] = 0;
        }
        if (target[1] + _C[q][1] == 0) {
          target[1] = _domainSizeY + 1;
        } else if (target[1] + _C[q][1] == _domainSizeY + 1) {
          target[1] = 0;
        }
        if (target[2] + _C[q][2] == 0) {
          target[2] = _domainSizeZ + 1;
        } else if (target[2] + _C[q][2] == _domainSizeZ + 1) {
          target[2] = 0;
        }
        const int periodicNb = target[0] + (_domainSizeX + 2) * (target[1] + (_domainSizeY + 2) * target[2]);
        pdf[19 * periodicNb + q] = pdf[index + q];
      }
    }
  }

  /** @brief refers to the LB method; computes density and velocity on pdf
   *  @param vel velocity
   *  @param density density
   *  @param pdf partial distribution function */
  void computeDensityAndVelocity(double* const vel, double& density, const double* const pdf) {
    vel[0] = -(pdf[1] + pdf[5] + pdf[8] + pdf[11] + pdf[15]);
    density = pdf[3] + pdf[7] + pdf[10] + pdf[13] + pdf[17];
    vel[1] = (pdf[4] + pdf[11] + pdf[12] + pdf[13] + pdf[18]) - (pdf[0] + pdf[5] + pdf[6] + pdf[7] + pdf[14]);
    vel[0] = density + vel[0];
    density = density + pdf[0] + pdf[1] + pdf[2] + pdf[4] + pdf[5] + pdf[6] + pdf[8] + pdf[9] + pdf[11] + pdf[12] + pdf[14] + pdf[15] + pdf[16] + pdf[18];
    vel[2] = (pdf[14] + pdf[15] + pdf[16] + pdf[17] + pdf[18]) - (pdf[0] + pdf[1] + pdf[2] + pdf[3] + pdf[4]);
    vel[0] = vel[0] / density;
    vel[1] = vel[1] / density;
    vel[2] = vel[2] / density;
  }

  /** takes care of communication across one face in one direction.
   *  @param pdf partial distribution function
   *  @param sendBuffer send buffer
   *  @param recvBuffer receive buffer
   *  @param nbFlagTo direction into which message is sent
   *  @param nbFlagFrom direction from which message is received
   *  @param startSend 3d coordinates that define the start of the data to be
   * sent to neighbouring process
   *  @param endSend 3d coordinates that define the end of the data to to be
   * sent to neighbouring process
   *  @param startRecv 3d coordinates that define the start of the data to be
   * received from neighbouring process
   *  @param endRecv 3d coordinates that define the end of the data to be
   * received from neighbouring process */
  void communicatePart(double* pdf, double* sendBuffer, double* recvBuffer, NbFlag nbFlagTo, NbFlag nbFlagFrom, tarch::la::Vector<3, int> startSend,
                       tarch::la::Vector<3, int> endSend, tarch::la::Vector<3, int> startRecv, tarch::la::Vector<3, int> endRecv) {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    // directions that point to LEFT/RIGHT,... -> same ordering as enums!
    const int directions[6][5] = {{1, 5, 8, 11, 15}, {3, 7, 10, 13, 17}, {4, 11, 12, 13, 18}, {0, 5, 6, 7, 14}, {0, 1, 2, 3, 4}, {14, 15, 16, 17, 18}};
    MPI_Request requests[2];
    MPI_Status status[2];
    tarch::la::Vector<2, int> plane;
    tarch::la::Vector<2, int> domainSize;
    // find out plane coordinates
    if (nbFlagTo == LEFT || nbFlagTo == RIGHT) {
      plane[0] = 1;
      plane[1] = 2;
      domainSize[0] = _domainSizeY;
      domainSize[1] = _domainSizeZ;
    } else if (nbFlagTo == FRONT || nbFlagTo == BACK) {
      plane[0] = 0;
      plane[1] = 2;
      domainSize[0] = _domainSizeX;
      domainSize[1] = _domainSizeZ;
    } else if (nbFlagTo == TOP || nbFlagTo == BOTTOM) {
      plane[0] = 0;
      plane[1] = 1;
      domainSize[0] = _domainSizeX;
      domainSize[1] = _domainSizeY;
    } else {
      std::cout << "ERROR LBCouetteSolver::communicatePart: d >2 or d < 0!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // extract data and write to send buffer
    tarch::la::Vector<3, int> coords(0);
    for (coords[2] = startSend[2]; coords[2] < endSend[2]; coords[2]++) {
      for (coords[1] = startSend[1]; coords[1] < endSend[1]; coords[1]++) {
        for (coords[0] = startSend[0]; coords[0] < endSend[0]; coords[0]++) {
          for (int q = 0; q < 5; q++) {
            sendBuffer[q + 5 * getParBuf(coords[plane[0]], coords[plane[1]], domainSize[0], domainSize[1])] =
                pdf[directions[nbFlagTo][q] + 19 * get(coords[0], coords[1], coords[2])];
          }
        }
      }
    }
    // send and receive data
    MPI_Irecv(recvBuffer, (domainSize[0] + 2) * (domainSize[1] + 2) * 5, MPI_DOUBLE, _parallelNeighbours[nbFlagFrom], 1000, coupling::indexing::IndexingService<3>::getInstance().getComm(), &requests[0]);
    MPI_Isend(sendBuffer, (domainSize[0] + 2) * (domainSize[1] + 2) * 5, MPI_DOUBLE, _parallelNeighbours[nbFlagTo], 1000, coupling::indexing::IndexingService<3>::getInstance().getComm(), &requests[1]);
    MPI_Waitall(2, requests, status);
    // write data back to pdf field
    if (_parallelNeighbours[nbFlagFrom] != MPI_PROC_NULL) {
      for (coords[2] = startRecv[2]; coords[2] < endRecv[2]; coords[2]++) {
        for (coords[1] = startRecv[1]; coords[1] < endRecv[1]; coords[1]++) {
          for (coords[0] = startRecv[0]; coords[0] < endRecv[0]; coords[0]++) {
            for (int q = 0; q < 5; q++) {
              if (_flag[get(coords[0], coords[1], coords[2])] == PARALLEL_BOUNDARY) {
                pdf[directions[nbFlagTo][q] + 19 * get(coords[0], coords[1], coords[2])] =
                    recvBuffer[q + 5 * getParBuf(coords[plane[0]], coords[plane[1]], domainSize[0], domainSize[1])];
              }
            }
          }
        }
      }
    }
#endif
  }

  /** @brief comunicates the boundary field data between the different processes
   */
  void communicate() {
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    // send from right to left
    communicatePart(_pdf1, _sendBufferX, _recvBufferX, LEFT, RIGHT, tarch::la::Vector<3, int>(1, 1, 1),
                    tarch::la::Vector<3, int>(2, _domainSizeY + 1, _domainSizeZ + 1), tarch::la::Vector<3, int>(_domainSizeX + 1, 1, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 1, _domainSizeZ + 1));
    // send from left to right
    communicatePart(_pdf1, _sendBufferX, _recvBufferX, RIGHT, LEFT, tarch::la::Vector<3, int>(_domainSizeX, 1, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 1, _domainSizeY + 1, _domainSizeZ + 1), tarch::la::Vector<3, int>(0, 1, 1),
                    tarch::la::Vector<3, int>(1, _domainSizeY + 1, _domainSizeZ + 1));
    // send from back to front
    communicatePart(_pdf1, _sendBufferY, _recvBufferY, FRONT, BACK, tarch::la::Vector<3, int>(0, 1, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, 2, _domainSizeZ + 1), tarch::la::Vector<3, int>(0, _domainSizeY + 1, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 2, _domainSizeZ + 1));
    // send from front to back
    communicatePart(_pdf1, _sendBufferY, _recvBufferY, BACK, FRONT, tarch::la::Vector<3, int>(0, _domainSizeY, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 1, _domainSizeZ + 1), tarch::la::Vector<3, int>(0, 0, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, 1, _domainSizeZ + 1));
    // send from top to bottom
    communicatePart(_pdf1, _sendBufferZ, _recvBufferZ, BOTTOM, TOP, tarch::la::Vector<3, int>(0, 0, 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 2, 2), tarch::la::Vector<3, int>(0, 0, _domainSizeZ + 1),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 2, _domainSizeZ + 2));
    // send from bottom to top
    communicatePart(_pdf1, _sendBufferZ, _recvBufferZ, TOP, BOTTOM, tarch::la::Vector<3, int>(0, 0, _domainSizeZ),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 2, _domainSizeZ + 1), tarch::la::Vector<3, int>(0, 0, 0),
                    tarch::la::Vector<3, int>(_domainSizeX + 2, _domainSizeY + 2, 1));
#endif
  }

  /** @brief relaxation frequency */
  const double _omega;
  /** @brief velocity of moving wall of Couette flow */
  tarch::la::Vector<3, double> _wallVelocity;
  /** @brief partical distribution function field */
  double* _pdf1{NULL};
  /** @brief partial distribution function field (stores the old time step)*/
  double* _pdf2{NULL};
  /** @brief lattice velocities*/
  const int _C[19][3]{{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0},
                      {1, 0, 0},   {-1, 1, 0},  {0, 1, 0},  {1, 1, 0},  {0, -1, 1}, {-1, 0, 1},  {0, 0, 1},  {1, 0, 1},  {0, 1, 1}};
  /** @brief lattice weights */
  const double _W[19]{1.0 / 36.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 3.0,
                      1.0 / 18.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0};
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVER_H_
