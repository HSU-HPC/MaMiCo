// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_
#define _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
namespace datastructures {
template <unsigned int dim> class CouplingCells;
} // namespace datastructures
} // namespace coupling

/**
 *	@brief provides access to the coupling cells. Base class for the class
 *coupling::datastructures::CouplingCellsWithLinkedCells
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::datastructures::CouplingCells {
public:
  /** Constructor: initialises the coupling cell
   * 	@param mdSolverInterface
   */
  CouplingCells(const coupling::IndexConversion<dim>& indexConversion);
  /** Destructor */
  ~CouplingCells();

  /** returns vector-of-pointers to coupling cells without access to linked
   * cells. We use this structure for data exchange between macroscopic and MD
   * solver.
   */
  const std::vector<coupling::datastructures::CouplingCell<dim>*>& getCouplingCells() const;

protected:
  /** holds pointers to all coupling cells with linked cells, but without
   * access to linked cells. This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
};

#include "CouplingCells.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_DATASTRUCTURES_COUPLINGCELLS_H_
