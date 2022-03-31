// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _DUMMY_SOLVER_INTERFACE_H_
#define _DUMMY_SOLVER_INTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

/** Macroscopic Solver Interface Implementation for Dummy Solver.
 *  The dummy solver always runs on rank 0 only.
 *  @author Rahul Arora
 */

class DummySolverInterface
    : public coupling::interface::MacroscopicSolverInterface<3> {
public:
  DummySolverInterface(
      tarch::la::Vector<3, unsigned int> globalNumberMacroscopicCells)
      : _globalNumberMacroscopicCells(globalNumberMacroscopicCells),
        _numberOfBoundaryCells(3) {}

  ~DummySolverInterface() {}

  bool receiveMacroscopicQuantityFromMDSolver(
      tarch::la::Vector<3, unsigned int> globalCellIndex) {
    bool isInside = true;
    for (int d = 0; d < 3; d++) {
      isInside =
          isInside && (globalCellIndex[d] > _numberOfBoundaryCells - 1) &&
          (globalCellIndex[d] <
           _globalNumberMacroscopicCells[d] + 2 - _numberOfBoundaryCells);
    }
    return isInside;
  }

  bool sendMacroscopicQuantityToMDSolver(
      tarch::la::Vector<3, unsigned int> globalCellIndex) {
    return (!receiveMacroscopicQuantityFromMDSolver(globalCellIndex));
  }

  std::vector<unsigned int>
  getRanks(tarch::la::Vector<3, unsigned int> globalCellIndex) {
    std::vector<unsigned int> ranks;
    ranks.push_back(0);
    return ranks;
  }

private:
  const tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
  const unsigned int _numberOfBoundaryCells;
};
#endif // _DUMMY_SOLVER_INTERFACE_H
