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
// each operation. this currently supports a maximum of <10M instances of
// MacroscopicCellService.
#define TAG_FROM_MD2MACRO 100000
#define TAG_FROM_MACRO2MD 500000

#include "tarch/la/Vector.h"
#include "tarch/utils/OstreamOperators.h"
#include <map>
#include <set>
#include <vector>

namespace coupling {
/** This is an enum for the macroscopic solver ID
 * @todo I have no idea what this is @someone */
enum MacroscopicSolverID {
  PEANO_LATTICEBOLTZMANN_ID = 0,       ///< test1
  TEST_LOCAL_MACROSCOPIC_SOLVER_ID = 1 ///< test2
};

// --------------------------- HELPER FUNCTIONS
// ----------------------------------------

/** @brief initialises the range for looping over Cartesian grid for general
 * case
 *  @param vec the input vector to be generalised
 *  @tparam dim refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @returns the given vector as 3d version */
template <unsigned int dim> tarch::la::Vector<3, unsigned int> initRange(tarch::la::Vector<dim, unsigned int> vec) {
  tarch::la::Vector<3, unsigned int> range(1);
  for (unsigned int d = 0; d < dim; d++) {
    range[d] = vec[d];
  }
  return range;
}

/** @brief reduces three-dimensional vector to vector of size dim (assumes that
 * dim<= 3).
 *  @param vec the input which will be reduced
 *  @tparam dim refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @return the given 3d vector as a vector for the simulation dimension (dim)
 */
template <unsigned int dim> tarch::la::Vector<dim, unsigned int> initDimVector(tarch::la::Vector<3, unsigned int> vec) {
  tarch::la::Vector<dim, unsigned int> smallVec;
  for (unsigned int d = 0; d < dim; d++) {
    smallVec[d] = vec[d];
  }
  return smallVec;
}

/** @briefcomputes a vectorwise division factor for the conversion of linearised
 * to vector indices
 *  @param numberCells the number of cells per direction in the spacial domain
 *  @tparam dim refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @returns the division factor for the simulation */
template <unsigned int dim> tarch::la::Vector<dim, unsigned int> initDivisionFactor(tarch::la::Vector<dim, unsigned int> numberCells) {
  tarch::la::Vector<dim, unsigned int> divFactor(1);
  for (unsigned int d = 1; d < dim; d++) {
    divFactor[d] = divFactor[d - 1] * (numberCells[d - 1]);
  }
  return divFactor;
}

/** This method is supposed to be only called with member divisionFactor
 * variables, cf. initDivisionFactor()-method.
 *  @brief converts linearised cell index to a vector cell index using
 * predefined division factors.
 *  @param cellIndex the index of the cell which shall be transformed
 *  @param divisionFactor the division factor for the corresponding spacial
 * setup of the domain
 *  @tparam dim refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @returns the dimensionized index for a given cell */
template <unsigned int dim>
tarch::la::Vector<dim, unsigned int> getVectorCellIndex(unsigned int cellIndex, const tarch::la::Vector<dim, unsigned int>& divisionFactor) {
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
