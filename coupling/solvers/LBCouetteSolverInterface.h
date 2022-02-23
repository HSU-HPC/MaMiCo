// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVERINTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
namespace solvers {
class LBCouetteSolverInterface;
}
} // namespace coupling

/** We only receive data from MD in the inner region, and we only send data for
 *  the outer region to MD. What "outer" means is specified by the arguments of
 * the interface; moreover, we do not send the ghost layer data from Couette to
 * MD. "inner" refers to "not outer" ;-) By default, the argument outerRegion is
 * set to 1: this yields that only the first non-ghost layer of macroscopic
 * cells shall be sent to MD, and all other inner macroscopic cells are received
 * from MD.
 *  @brief interface for the LBCouetteSolver
 *  @author Philipp Neumann  */
class coupling::solvers::LBCouetteSolverInterface
    : public coupling::interface::MacroscopicSolverInterface<3> {
public:
  /** @brief a simple constructor
   *  @param avgNumberLBCells the average number of cells per process of the
   * lattice Boltzmann solver (dimensioned)
   *  @param numberProcesses the total number of mpi processes on which the
   * solver is parallelised
   *  @param offsetMDDomain offset (measured in cell units) of the MD domain
   * (excluding any (LB/MD) ghost layers)
   *  @param globalNumberMacroscopicCells the total number of macroscopic cells
   *  @param outerRegion defines, how many cell layers will be sent to the macro
   * solver */
  LBCouetteSolverInterface(
      tarch::la::Vector<3, unsigned int> avgNumberLBCells,
      tarch::la::Vector<3, unsigned int> numberProcesses,
      tarch::la::Vector<3, unsigned int> offsetMDDomain,
      tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells,
      unsigned int outerRegion = 1)
      : _avgNumberLBCells(avgNumberLBCells), _numberProcesses(numberProcesses),
        _offsetMDDomain(offsetMDDomain), _outerRegion(outerRegion),
        _globalNumberMacroscopicCells(globalNumberMacroscopicCells) {}
  ~LBCouetteSolverInterface() {}

  /** with this function one can check, which data needs so be send from micro
   * to macro solver for the correct Couette scenario setup, all (inner) cells
   * need to be received
   *  @brief checks for a given macroscopic cell, if it needs to be received
   * (true) ro not (false)
   *  @param globalCellIndex global dimensioned cell index to check for
   *  @returns a bool, which indicates if the cell will be received*/
  bool receiveMacroscopicQuantityFromMDSolver(
      tarch::la::Vector<3, unsigned int> globalCellIndex) {
    bool recv = true;
    for (unsigned int d = 0; d < 3; d++) {
      recv = recv && (globalCellIndex[d] > _outerRegion) &&
             (globalCellIndex[d] <
              _globalNumberMacroscopicCells[d] + 1 - _outerRegion);
    }
    return recv;
  }

  /** send all macroscopic cell data within a boundary strip to MD. Only send
   * data that are not in the ghost layer and not part of the inner region.
   *  @brief checks for a given cell if it needs to be send (true) or not
   * (false)
   *  @param globalCellIndex global dimensioned cell index to check for
   *  @returns a bool, which indicates if the cell will be send */
  bool sendMacroscopicQuantityToMDSolver(
      tarch::la::Vector<3, unsigned int> globalCellIndex) {
    bool outer = false;
    for (unsigned int d = 0; d < 3; d++) {
      outer = outer || (globalCellIndex[d] < 1) ||
              (globalCellIndex[d] > _globalNumberMacroscopicCells[d]);
    }
    return (!outer) &&
           (!receiveMacroscopicQuantityFromMDSolver(globalCellIndex));
  }

  /** @brief returns for a given macroscopic cell index, which rank holds the
   * correct data
   *  @oaram globalCellIndex global dimensioned cell index to check for
   *  @returns a vector containing all correct ranks  */
  virtual std::vector<unsigned int>
  getRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    std::vector<unsigned int> ranks;
    // determine global index of cell in LB simulation
    tarch::la::Vector<3, unsigned int> globalLBCellIndex(globalCellIndex +
                                                         _offsetMDDomain);
    // modify global LB cell index due to ghost layer
    for (int d = 0; d < 3; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "LB cell index for global cell index " << globalCellIndex
                << ": " << globalLBCellIndex << std::endl;
#endif
      if (globalLBCellIndex[d] > 0) {
        globalLBCellIndex[d]--;
      }
    }
    // loop over all neighbouring cells within a one-cell surrounding and detect
    // the respective ranks. IMPROVE: This currently only allows for simulations
    // with MD located inside the domain (no simulation across boundary)
    for (int z = -1; z < 2; z++) {
      for (int y = -1; y < 2; y++) {
        for (int x = -1; x < 2; x++) {
          // neighbour cell index
          const tarch::la::Vector<3, unsigned int> nbIndex(
              globalLBCellIndex[0] + x, globalLBCellIndex[1] + y,
              globalLBCellIndex[2] + z);
          // coordinates of process of neighbour celll
          const tarch::la::Vector<3, unsigned int> processCoordinates(
              (nbIndex[0] / _avgNumberLBCells[0]),
              (nbIndex[1] / _avgNumberLBCells[1]),
              (nbIndex[2] / _avgNumberLBCells[2]));
          // corresponding rank
          const unsigned int rank =
              processCoordinates[0] +
              _numberProcesses[0] *
                  (processCoordinates[1] +
                   processCoordinates[2] * _numberProcesses[1]);

          // if this rank is not part of the vector, push it back
          bool found = false;
          for (unsigned int i = 0; i < ranks.size(); i++) {
            found = found || (rank == ranks[i]);
          }
          if (!found) {
            ranks.push_back(rank);
          }
        }
      }
    }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Ranks for cell " << globalCellIndex << ":";
    for (unsigned int i = 0; i < ranks.size(); i++) {
      std::cout << " " << ranks[i];
    }
    std::cout << std::endl;
#endif
    return ranks;
  }

  /** provides the ranks of the Couette solver that send valid LB flow data to
   * the coupling tool. Since we may have multiple copies of a cell due to the
   * ghost layers of the LBCouetteSolver, those ghost cell copies do not contain
   * valid information. We therefore need to come up with a special
   * implementation for getSourceRanks(). This implementation only returns one
   * rank per cell, that is the rank which holds the non-ghost cell copy.
   *  @brief returns for a given macroscopic cell index, which source rank holds
   * the correct data
   *  @param globalCellIndex global dimensioned cell index to check for
   *  @returns the vector of the correct rank  */
  virtual std::vector<unsigned int>
  getSourceRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    // determine global index of cell in LB simulation
    tarch::la::Vector<3, unsigned int> globalLBCellIndex(globalCellIndex +
                                                         _offsetMDDomain);
    // modify global LB cell index due to ghost layer
    for (int d = 0; d < 3; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "LB cell index for global cell index " << globalCellIndex
                << ": " << globalLBCellIndex << std::endl;
#endif
      if (globalLBCellIndex[d] > 0) {
        globalLBCellIndex[d]--;
      }
    }
    // determine process coordinates and respective rank
    const tarch::la::Vector<3, unsigned int> processCoordinates(
        globalLBCellIndex[0] / _avgNumberLBCells[0],
        globalLBCellIndex[1] / _avgNumberLBCells[1],
        globalLBCellIndex[2] / _avgNumberLBCells[2]);
    const unsigned int rank =
        processCoordinates[0] +
        _numberProcesses[0] * (processCoordinates[1] +
                               processCoordinates[2] * _numberProcesses[1]);
    std::vector<unsigned int> ranks;
    ranks.push_back(rank);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Source rank for cell " << globalCellIndex << ": " << ranks[0]
              << std::endl;
#endif
    return ranks;
  }

private:
  /** @brief avg. number of LB cells per LB process (must be same for Interface
   * and LBCouetteSolver) */
  const tarch::la::Vector<3, unsigned int> _avgNumberLBCells;
  /** @brief number of processes used by LB solver */
  const tarch::la::Vector<3, unsigned int> _numberProcesses;
  /** @brief offset of MD domain (excl. any ghost layers on MD or LB side) */
  const tarch::la::Vector<3, unsigned int> _offsetMDDomain;
  /** @brief defines an offset of cells which is considered to be the outer
   * region */
  const unsigned int _outerRegion;
  /** @brief global number of macroscopic cells */
  const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVERINTERFACE_H_
