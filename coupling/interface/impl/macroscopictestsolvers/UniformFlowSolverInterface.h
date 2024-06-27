// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_UNIFORMFLOWSOLVERINTERFACE_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_UNIFORMFLOWSOLVERINTERFACE_H_

#include "coupling/interface/impl/macroscopictestsolvers/TestMacroscopicSolverInterface.h"

namespace coupling {
namespace interface {
template <unsigned int dim> class UniformFlowSolverInterface;
}
} // namespace coupling

/** accelerates the fluid along x-direction with values provided everywhere. No values are sent back from MD to this (kind of) flow solver.
 *  The cells of this solver are assumed to be located on rank 0.
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::interface::UniformFlowSolverInterface : public coupling::interface::TestMacroscopicSolverInterface<dim> {
public:
  UniformFlowSolverInterface(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells, double massPerCouplingCell,
                             tarch::la::Vector<dim, double> momentumPerCouplingCell)
      : TestMacroscopicSolverInterface<dim>(), _totalNumberCouplingCells(initTotalNumberCouplingCells(globalNumberCouplingCells)),
        _massPerCouplingCell(massPerCouplingCell), _momentumPerCouplingCell(momentumPerCouplingCell) {}
  virtual ~UniformFlowSolverInterface() {}

  /** no values are received. */
  virtual bool receiveMacroscopicQuantityFromMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return false; }

  /** In every grid cell, values are sent.
   */
  virtual bool sendMacroscopicQuantityToMDSolver(tarch::la::Vector<dim, unsigned int> globalCellIndex) { return true; }

  /** all cells are located on rank 0.
   */
  virtual std::vector<unsigned int> getRanks(tarch::la::Vector<dim, unsigned int> globalCellIndex) {
    std::vector<unsigned int> result;
    result.push_back(0);
    return result;
  }

  /** add every coupling cell with mass and momentum to the send buffer */
  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Sending() {
    std::vector<coupling::datastructures::CouplingCell<dim>*> result;
    for (unsigned int i = 0; i < _totalNumberCouplingCells; i++) {
      result.push_back(new coupling::datastructures::CouplingCell<dim>());
      if (result[i] == NULL) {
        std::cout << "ERROR coupling::interface::UniformFlowSolverInterface: result[i]==NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
      result[i]->setMicroscopicMass(_massPerCouplingCell);
      result[i]->setMicroscopicMomentum(_momentumPerCouplingCell);
    }
    return result;
  }
  /** add every index to the send buffer */
  virtual unsigned int* getCouplingCellIndices4Sending() {
    unsigned int* indices = new unsigned int[_totalNumberCouplingCells];
    if (indices == NULL) {
      std::cout << "ERROR coupling::interface::UniformFlowSolverInterface: indices==NULL!" << std::endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < _totalNumberCouplingCells; i++) {
      indices[i] = i;
    }
    return indices;
  }

  /** nothing is received */
  virtual std::vector<coupling::datastructures::CouplingCell<dim>*> getCouplingCells4Receiving() {
    std::vector<coupling::datastructures::CouplingCell<dim>*> result;
    return result;
  }
  virtual unsigned int* getCouplingCellIndices4Receiving() { return NULL; }

private:
  /** total number of coupling cells incl. ghost layers */
  const unsigned int _totalNumberCouplingCells;
  /** mass and momentum that is assumed to be placed in each coupling cell on this solver's side. */
  const double _massPerCouplingCell;
  const tarch::la::Vector<dim, double> _momentumPerCouplingCell;

  /** returns the total number of coupling cells incl. ghost layers */
  unsigned int initTotalNumberCouplingCells(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells) const {
    // total number cells = global number cells + ghost layers
    tarch::la::Vector<dim, unsigned int> totalNumberCells = globalNumberCouplingCells + tarch::la::Vector<dim, unsigned int>(2);
    unsigned int num = totalNumberCells[0];
    for (unsigned int d = 1; d < dim; d++) {
      num = num * totalNumberCells[d];
    }
    return num;
  }
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_UNIFORMFLOWSOLVERINTERFACE_H_
