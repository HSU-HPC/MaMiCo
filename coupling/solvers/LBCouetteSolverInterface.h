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

/** interface to couette solver. We only receive data from MD in the inner
 * region, and we only send data for the outer region to MD. What "outer" means
 * is specified by the arguments of the interface;
 *  moreover, we do not send the ghost layer data from Couette to MD.
 *  "inner" refers to "not outer" ;-)
 *  By default, the argument outerRegion is set to 1: this yields that only the
 * first non-ghost layer of macroscopic cells shall be sent to MD, and all other
 * inner macroscopic cells are received from MD.
 *  @author Philipp Neumann
 */
class coupling::solvers::LBCouetteSolverInterface : public coupling::interface::MacroscopicSolverInterface<3> {
public:
  LBCouetteSolverInterface(tarch::la::Vector<3, unsigned int> avgNumberLBCells, // avg number of LB cells per process (obtained from
                                                                                // LBCouetteSolver)
                           tarch::la::Vector<3, unsigned int> numberProcesses,  // total number of processes used by the LB Couette
                                                                                // solver (obtained from LBCouetteSolver)
                           tarch::la::Vector<3, unsigned int> offsetMDDomain,   // offset (measured in cell units) of the MD domain
                                                                                // (excluding any (LB/MD) ghost layers)
                           tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells,
                           unsigned int outerRegion = 1)
      : // see CouetteSolverInterface
        _avgNumberLBCells(avgNumberLBCells), _numberProcesses(numberProcesses), _offsetMDDomain(offsetMDDomain), _outerRegion(outerRegion),
        _globalNumberMacroscopicCells(globalNumberMacroscopicCells) {}
  ~LBCouetteSolverInterface() {}

  /** receive all (inner) cells */
  bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    bool recv = true;
    for (unsigned int d = 0; d < 3; d++) {
      recv = recv && (globalCellIndex[d] > _outerRegion) && (globalCellIndex[d] < _globalNumberMacroscopicCells[d] + 1 - _outerRegion);
    }
    return recv;
  }

  /** send all macroscopic cell data within a boundary strip to MD. Only send
   * data that are not in the ghost layer and not part of the inner region. */
  bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    bool outer = false;
    for (unsigned int d = 0; d < 3; d++) {
      outer = outer || (globalCellIndex[d] < 1) || (globalCellIndex[d] > _globalNumberMacroscopicCells[d]);
    }
    return (!outer) && (!receiveMacroscopicQuantityFromMDSolver(globalCellIndex));
  }

  /** returns all ranks that the cell at globalCellIndex is associated to. We
   * use this definition to receive data (cf. method getTargetRanks() in the
   * MacroscopicSolverInterface)
   *  since we may also require valid data for macroscopic cells in the ghost
   * layers of the LBCouetteSolver.
   */
  virtual std::vector<unsigned int> getRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    std::vector<unsigned int> ranks;

    // determine global index of cell in LB simulation
    tarch::la::Vector<3, unsigned int> globalLBCellIndex(globalCellIndex + _offsetMDDomain);
    // modify global LB cell index due to ghost layer
    for (int d = 0; d < 3; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "LB cell index for global cell index " << globalCellIndex << ": " << globalLBCellIndex << std::endl;
#endif
      if (globalLBCellIndex[d] > 0) {
        globalLBCellIndex[d]--;
      }
    }

    // loop over all neighbouring cells within a one-cell surrounding and detect
    // the respective ranks. IMPROVE: This currently only allows for simulations
    // with MD located inside the domain
    // (no simulation across boundary)
    for (int z = -1; z < 2; z++) {
      for (int y = -1; y < 2; y++) {
        for (int x = -1; x < 2; x++) {
          // neighbour cell index
          const tarch::la::Vector<3, unsigned int> nbIndex(globalLBCellIndex[0] + x, globalLBCellIndex[1] + y, globalLBCellIndex[2] + z);
          // coordinates of process of neighbour celll
          const tarch::la::Vector<3, unsigned int> processCoordinates((nbIndex[0] / _avgNumberLBCells[0]), (nbIndex[1] / _avgNumberLBCells[1]),
                                                                      (nbIndex[2] / _avgNumberLBCells[2]));
          // corresponding rank
          const unsigned int rank = processCoordinates[0] + _numberProcesses[0] * (processCoordinates[1] + processCoordinates[2] * _numberProcesses[1]);

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
   * ghost layers of
   *  the LBCouetteSolver, those ghost cell copies do not contain valid
   * information. We therefore need to come up with a special implementation for
   * getSourceRanks().
   *  This implementation only returns one rank per cell, that is the rank which
   * holds the non-ghost cell copy. */
  virtual std::vector<unsigned int> getSourceRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    // determine global index of cell in LB simulation
    tarch::la::Vector<3, unsigned int> globalLBCellIndex(globalCellIndex + _offsetMDDomain);
    // modify global LB cell index due to ghost layer
    for (int d = 0; d < 3; d++) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
      std::cout << "LB cell index for global cell index " << globalCellIndex << ": " << globalLBCellIndex << std::endl;
#endif
      if (globalLBCellIndex[d] > 0) {
        globalLBCellIndex[d]--;
      }
    }
    // determine process coordinates and respective rank
    const tarch::la::Vector<3, unsigned int> processCoordinates(globalLBCellIndex[0] / _avgNumberLBCells[0], globalLBCellIndex[1] / _avgNumberLBCells[1],
                                                                globalLBCellIndex[2] / _avgNumberLBCells[2]);
    const unsigned int rank = processCoordinates[0] + _numberProcesses[0] * (processCoordinates[1] + processCoordinates[2] * _numberProcesses[1]);
    std::vector<unsigned int> ranks;
    ranks.push_back(rank);

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "Source rank for cell " << globalCellIndex << ": " << ranks[0] << std::endl;
#endif
    return ranks;
  }

private:
  const tarch::la::Vector<3, unsigned int> _avgNumberLBCells;             // avg. number of LB cells per LB process (must be same
                                                                          // for Interface and LBCouetteSolver)
  const tarch::la::Vector<3, unsigned int> _numberProcesses;              // number of processes used by LB solver
  const tarch::la::Vector<3, unsigned int> _offsetMDDomain;               // offset of MD domain (excl. any ghost layers on MD or
                                                                          // LB side)
  const unsigned int _outerRegion;                                        // defines an offset of cells which is
                                                                          // considered to be the outer region
  const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells; // global number of macroscopic cells
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_LBCOUETTESOLVERINTERFACE_H_
