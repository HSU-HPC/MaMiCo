// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#ifndef _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_

#include "coupling/interface/MacroscopicSolverInterface.h"

namespace coupling {
  namespace solvers {
    template<unsigned int dim>
    class CouetteSolverInterface;
  }
}

/** We only receive data from MD in the inner region, and we only send data for the
 *  outer region to MD. What "outer" means is specified by the arguments of the interface;
 *  moreover, we do not send the ghost layer data from Couette to MD.
 *  "inner" refers to "not outer" ;-)
 *  By default, the argument outerRegion is set to 1: this yields that only the
 *  first non-ghost layer of macroscopic cells shall be sent to MD, and all other
 *  inner macroscopic cells are received from MD.
 *  The couette solver is expected to only run on rank 0.
 *  @brief interface to couette solver
 *  @author Philipp Neumann
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2, or 3 */
template<unsigned int dim>
class coupling::solvers::CouetteSolverInterface: public coupling::interface::MacroscopicSolverInterface<dim> {
  public:
    /** @brief a simple constructor
     *  @param globalNumberMacroscopicCells the total number of macroscopic cells
     *  @param outerRegion defines, how many cell layers will be sent to the macro solver */
    CouetteSolverInterface(tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,unsigned int outerRegion=1):
    coupling::interface::MacroscopicSolverInterface<dim>(),
    _outerRegion(outerRegion),_globalNumberMacroscopicCells(globalNumberMacroscopicCells){
    }

    /** @brief a dummy destructor */
    virtual ~CouetteSolverInterface(){}

    /** with this function one can check, which data needs so be send from micro to macro solver
     *  for the correct Couette scenario setup, all (inner) cells need to be received
     *  @brief checks for a given macroscopic cell, if it needs to be received or not
     *  @param globalCellIndex global dimensioned cell index to check for
     *  @returns a bool, which indicates if a cell will be received (true) or not (false) */
    virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim,unsigned int> globalCellIndex){
      bool recv=true;
      for (unsigned int d = 0; d < dim; d++){ recv = recv && (globalCellIndex[d]>_outerRegion) && (globalCellIndex[d]<_globalNumberMacroscopicCells[d]+1-_outerRegion); }
      return recv;
    }

    /** send all macroscopic cell data within a boundary strip to MD. Only send data
     *  that are not in the ghost layer and not part of the inner region.
     *  @brief checks for a given cell if it needs to be send or not
     *  @param globalCellIndex global dimensioned cell index to check for
     *  @returns a bool, which indicates if a cell will be send (true) or not (false)  */
    virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim,unsigned int> globalCellIndex){
      bool outer=false;
      for (unsigned int d = 0; d < dim; d++){ outer=outer || (globalCellIndex[d]<1) || (globalCellIndex[d]>_globalNumberMacroscopicCells[d]); }
      return (!outer) && (!receiveMacroscopicQuantityFromMDSolver(globalCellIndex));
    }

    /** @brief calculates for a macroscopic cell index, which rank holds it
     *  @param globalCellIndex global dimensioned cell index to check for
     *  @returns returns the rank on which the data is located */
    virtual std::vector<unsigned int> getRanks(tarch::la::Vector<dim,unsigned int> globalCellIndex){
      std::vector<unsigned int> ranks;
      ranks.push_back(0);
      return ranks;
    }

  private:
    /** @brief defines an offset of cells which is considered to be the outer region */
    const unsigned int _outerRegion;
    /** @brief global number of macroscopic cells */
    const tarch::la::Vector<dim,unsigned int> _globalNumberMacroscopicCells;
};

#endif // _MOLECULARDYNAMICS_COUPLING_SOLVERS_COUETTESOLVERINTERFACE_H_
