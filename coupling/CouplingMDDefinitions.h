// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_COUPLINGMDDEFINITIONS_H_
#define _MOLECULARDYNAMICS_COUPLING_COUPLINGMDDEFINITIONS_H_

#define COUPLING_MD_NO 0
#define COUPLING_MD_YES 1

#if defined(MDCoupledDebug)
#define COUPLING_MD_DEBUG COUPLING_MD_YES
#else
#define COUPLING_MD_DEBUG COUPLING_MD_NO
#endif

#if defined(MDCoupledError)
#define COUPLING_MD_ERROR COUPLING_MD_YES
#else
#define COUPLING_MD_ERROR COUPLING_MD_NO
#endif

#if defined(MDCoupledParallel)
#define COUPLING_MD_PARALLEL COUPLING_MD_YES
#else
#define COUPLING_MD_PARALLEL COUPLING_MD_NO
#endif
// message tags for MPI parallelisation. For each MacroscopicCellService, we add
// the ID of the cell service to the tags from below to yield a unique tag for
// each operation.
// this currently supports a maximum of <10M instances of
// MacroscopicCellService.
#define TAG_FROM_MD2MACRO 100000
#define TAG_FROM_MACRO2MD 500000

#if defined(MDCoupled)
#define COUPLING_MD_WITH_TEST_SOLVER COUPLING_MD_YES
#else
#define COUPLING_MD_WITH_TEST_SOLVER COUPLING_MD_NO
#endif

#include "tarch/la/Vector.h"

namespace coupling {
enum MacroscopicSolverID { PEANO_LATTICEBOLTZMANN_ID = 0, TEST_LOCAL_MACROSCOPIC_SOLVER_ID = 1 };

// --------------------------- HELPER FUNCTIONS
// ----------------------------------------

/** inits the range for looping over Cartesian grid for general case (1D/2D/3D
 * support).*/
template <unsigned int dim> tarch::la::Vector<3, unsigned int> initRange(tarch::la::Vector<dim, unsigned int> vec) {
  tarch::la::Vector<3, unsigned int> range(1);
  for (unsigned int d = 0; d < dim; d++) {
    range[d] = vec[d];
  }
  return range;
}

/** reduces three-dimensional vector to vector of size dim (assumes that dim<=
 * 3). */
template <unsigned int dim> tarch::la::Vector<dim, unsigned int> initDimVector(tarch::la::Vector<3, unsigned int> vec) {
  tarch::la::Vector<dim, unsigned int> smallVec;
  for (unsigned int d = 0; d < dim; d++) {
    smallVec[d] = vec[d];
  }
  return smallVec;
}

/** computes a vectorwise division factor for the conversion of linearised to
 * vector indices */
template <unsigned int dim> tarch::la::Vector<dim, unsigned int> initDivisionFactor(tarch::la::Vector<dim, unsigned int> numberCells) {
  tarch::la::Vector<dim, unsigned int> divFactor(1);
  for (unsigned int d = 1; d < dim; d++) {
    divFactor[d] = divFactor[d - 1] * (numberCells[d - 1]);
  }
  return divFactor;
}

/** converts linearised cell index to a vector cell index using predefined
 * division factors.
 *  This method is supposed to be only called with member divisionFactor
 * variables, cf. initDivisionFactor()-method.
 */
template <unsigned int dim>
tarch::la::Vector<dim, unsigned int> getVectorCellIndex(unsigned int cellIndex, const tarch::la::Vector<dim, unsigned int> &divisionFactor) {
  tarch::la::Vector<dim, unsigned int> myVector(0);
  unsigned int help = cellIndex;
  for (int d = dim - 1; d > 0; d--) {
    myVector[d] = help / divisionFactor[d];
    help = help - myVector[d] * divisionFactor[d];
  }
  myVector[0] = help;

  return myVector;
}

} // namespace coupling
#endif // _MOLECULARDYNAMICS_COUPLING_COUPLINGMDDEFINITIONS_H_
