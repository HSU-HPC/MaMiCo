// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace datastructures {
template <class CellIndexT, unsigned int dim> class CellContainer;
} // namespace datastructures
} // namespace coupling

/**
 *	@brief provides access to the coupling cells. Base class for the class
 *coupling::datastructures::LinkedCellContainer
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class CellIndexT, unsigned int dim> class coupling::datastructures::CellContainer {

public:
  /** returns a pointer to the coupling cell without access to linked cells. */
  const coupling::datastructures::CouplingCell<dim>* operator[](CellIndexT index) const;

  /** adds a new coupling cell to the datastructure at the next index (will only work if the data structure is not yet full)*/
  void operator<<(coupling::datastructures::CouplingCell<dim>* couplingCell);

protected:
  /** holds pointers to all coupling cells with linked cells, but without
   * access to linked cells. This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
};

#include "CellContainer.cpph"
