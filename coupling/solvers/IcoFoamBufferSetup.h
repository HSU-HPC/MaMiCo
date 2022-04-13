// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_ICOFOAMBUFFERSETUP_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_ICOFOAMBUFFERSETUP_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
namespace solvers {
template <unsigned int dim> class FoamSolverInterface;
}
} // namespace coupling

/** interface to the Foam solver. It's almost the same as the CouetteSolverInterface. That means, data from MD is received in the inner region plus one cell
 * layer from the outer region. Plus one is necessary cause the Foam solver needs the boundary condition on the boundary. So it is interpolated between the two
 * cells beside. The Foam solver sends data in tha outerRegion region, as defined by the parameters of the interface. . By default, the argument outerRegion is
 * set to 1: this yields that only the first non-ghost layer of macroscopic cells shall be sent to MD, and all other inner macroscopic cells are received from
 * MD. The Foam solver is expected to only run on rank 0.
 *  @author Philipp Neumann & Helene Wittenberg
 */
template <unsigned int dim> class coupling::solvers::FoamSolverInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  FoamSolverInterface(tarch::la::Vector<dim, unsigned int> globalNumberMacroscopicCells, unsigned int outerRegion = 1)
      : coupling::interface::MacroscopicSolverInterface<dim>(), _outerRegion(outerRegion), _globalNumberMacroscopicCells(globalNumberMacroscopicCells) {
  } virtual ~FoamSolverInterface() {
  }

  /** receive all (inner) cells */
  virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    bool recv = true;
    // for (unsigned int d = 0; d < dim; d++){ recv = recv && (globalCellIndex[d]>_outerRegion) &&
    // (globalCellIndex[d]<_globalNumberMacroscopicCells[d]+1-_outerRegion); }
    for (unsigned int d = 0; d < dim; d++) {
      recv = recv && (globalCellIndex[d] > _outerRegion - 1) && (globalCellIndex[d] < _globalNumberMacroscopicCells[d] + 2 - _outerRegion);
    }
    // for (unsigned int d = 0; d < dim; d++){ recv = recv && (globalCellIndex[d]>0) && (globalCellIndex[d]<_globalNumberMacroscopicCells[d]+1); }
    return recv;
  }

  /** send all macroscopic cell data within a boundary strip to MD. Only send data that are not in the ghost layer and not part of the inner region. */
  virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    bool outer = false;
    bool recv = true;
    for (unsigned int d = 0; d < dim; d++) {
      recv = recv && (globalCellIndex[d] > _outerRegion) && (globalCellIndex[d] < _globalNumberMacroscopicCells[d] + 1 - _outerRegion);
    }
    for (unsigned int d = 0; d < dim; d++) {
      outer = outer || (globalCellIndex[d] < 1) || (globalCellIndex[d] > _globalNumberMacroscopicCells[d]);
    }
    return (!outer) && (!recv);
  }

  virtual std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    std::vector<unsigned int> ranks;
    ranks.push_back(0);
    return ranks;
  }

private:
  const unsigned int _outerRegion;                                          // defines an offset of cells which is considered to be the outer region
  const tarch::la::Vector<dim, unsigned int> _globalNumberMacroscopicCells; // global number of macroscopic cells
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_ICOFOAMBUFFERSETUP_H_
