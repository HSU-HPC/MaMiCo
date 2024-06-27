// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
namespace solvers {
template <unsigned int dim> class CouetteSolverInterface;
}
} // namespace coupling

/** We only receive data from MD in the inner region, and we only send data for
 * the outer region to MD. What "outer" means is specified by the arguments of
 * the interface; moreover, we do not send the ghost layer data from Couette to
 * MD. "inner" refers to "not outer" ;-) By default, the argument outerRegion is
 * set to 1: this yields that only the first non-ghost layer of coupling
 * cells shall be sent to MD, and all other inner coupling cells are received
 * from MD. The couette solver is expected to only run on rank 0.
 *  @brief interface to couette solver
 *  @author Philipp Neumann
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3 */
template <unsigned int dim> class coupling::solvers::CouetteSolverInterface : public coupling::interface::MacroscopicSolverInterface<dim> {
public:
  /** @brief a simple constructor
   *  @param globalNumberCouplingCells the total number of coupling cells
   *  @param outerRegion defines, how many cell layers will be sent to the MD
   * solver */
  CouetteSolverInterface(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells, unsigned int outerRegion = 1)
      : coupling::interface::MacroscopicSolverInterface<dim>(), _outerRegion(outerRegion), _globalNumberCouplingCells(globalNumberCouplingCells) {}

  /** @brief a dummy destructor */
  virtual ~CouetteSolverInterface() {}

  unsigned int getOuterRegion() override { return _outerRegion; }

  /** @brief calculates for a coupling cell index, which rank holds it
   *  @param globalCellIndex global dimensioned cell index to check for
   *  @returns returns the rank on which the data is located */
  std::vector<unsigned int> getRanks(I01 idx) override {
    std::vector<unsigned int> ranks;
    ranks.push_back(0);
    return ranks;
  }

private:
  /** @brief defines an offset of cells which is considered to be the outer
   * region */
  const unsigned int _outerRegion;
  /** @brief global number of coupling cells */
  const tarch::la::Vector<dim, unsigned int> _globalNumberCouplingCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_
