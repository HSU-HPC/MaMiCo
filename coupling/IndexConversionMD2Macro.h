// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

// #define DEBUG_ICM2M
// #define DEBUG_ICM2M_VERBOSE

#include "IndexConversion.h"
#include "interface/MacroscopicSolverInterface.h"
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
#include <mpi.h>
#endif

namespace coupling {
template <unsigned int dim> class IndexConversionMD2Macro;
}

/** Wrapper class for coupling::IndexConversion for cases in which you need to
 * know boundaries of the MD-To-Macro (M2M) domain (that is the domain of cells
 * that are transfered from MD to continuum) Note that currently a rework of
 * indexing in MaMiCo takes place, making this obsolete in the future.
 *
 *  @author Felix Maurer
 */
template <unsigned int dim> class coupling::IndexConversionMD2Macro {
public:
  /** @brief a simple constructor
   *  @param indexConversion the index conversion class
   *  @param macroscopicSolverInterface the interface to the macroscopic solver
   * in application
   *  @param comm the mpi communicator to use
   *  @param lowestRankInComm the lowest rank in the mpi communicator
   */
  IndexConversionMD2Macro(const coupling::IndexConversion<dim>* indexConversion,
                          coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                          ,
                          const MPI_Comm comm, const int lowestRankInComm = 0
#endif
                          /** @todo @felix adapt the params comm and
                             lowestRankInComm for case multimd */
                          )
      : _ic(indexConversion), _msi(macroscopicSolverInterface), _lowerBoundaryAllRanks(nullptr), _upperBoundaryAllRanks(nullptr),
        _lowerBoundaryThisRank(nullptr), _upperBoundaryThisRank(nullptr)
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
        ,
        _comm(comm), _lowestRank((int)lowestRankInComm), _myRank(_ic != nullptr ? _ic->getThisRank() : -1)
#endif
  {
    if (_ic == nullptr)
      throw std::runtime_error("IndexConversionMD2Macro: Constructor called "
                               "with nullptr as base IndexConversion.");
    if (_msi == nullptr)
      throw std::runtime_error("IndexConversionMD2Macro: Constructor called "
                               "with nullptr as MacroscopicSolverInterface.");
#ifdef DEBUG_ICM2M
    std::cout << "ICM2M: Created new instance at location " << this << " using IC at " << _ic << " and MSI at " << _msi << std::endl;
#endif
  }

  /** @brief a simple destructor*/
  ~IndexConversionMD2Macro() {
    delete _lowerBoundaryAllRanks;
    delete _upperBoundaryAllRanks;
    delete _lowerBoundaryThisRank;
    delete _upperBoundaryThisRank;

#ifdef DEBUG_ICM2M
    std::cout << "ICM2M: Deconstructed." << std::endl;
#endif
  }

  /** Mimics initMD2MacroDomain's other implementation (see below) where you
   * already have the MD-To-Macro domain at hand. It then (trivially)
   * initializes Md-To-Macro and Outer cells/indices.
   *
   *  This function currently only has uses in cases where you do not have
   * access to all cells of the entire domain, which is only the case for post
   * multi-instance filtering as of now. In this case, the "outer" region is
   * empty.
   *
   *  @brief This method assumes the given cells to be the entire m2m-domain
   * (see below).
   *  @param m2mDomainCells Vector containing all MD-To-Macro cells.
   *  @param m2mGlobalCellIndices Vector containing all MD-To-Macro indices.
   *  @param outerCells Vector containing all outer cells.
   *  @param outerIndices Vector containing all outer indices.
   */
  void initMD2MacroDomain(std::vector<coupling::datastructures::CouplingCell<dim>*> m2mDomainCells,
                          std::vector<tarch::la::Vector<dim, unsigned int>> m2mGlobalCellIndices,
                          std::vector<coupling::datastructures::CouplingCell<dim>*>& outerCells,
                          std::vector<tarch::la::Vector<dim, unsigned int>>& outerIndices);

  /** This alternative chooses a subspace of the cell (and index) input based on
   * what will be transfered to the macro solver. This subspace is referred to
   * as "md2Macro-domain" (or sometimes (m2m-domain). All other cells will be
   * placed into the "outer" domain. Note that all indexing here is in terms of
   * MaMiCo-Indexing nonetheless, i.e. m2mDomainCells don't start counting at
   * the start of the MD-To-Macro domain.
   *
   *  @param inputCells Vector containg all cells in MD domain.
   *  @param m2mDomainCells Output vector where cells located within the
   * MD-To-Macro cells shall be placed.
   *  @param m2mIndices Output vector where indices of those MD-To-Macro cells
   * shall be placed.
   *  @param m2mDomainCells Output vector where cells located NOT within the
   * MD-To-Macro cells shall be placed.
   *  @param m2mIndices Output vector where indices of those outer cells shall
   * be placed.
   */
  void initMD2MacroDomain(std::vector<coupling::datastructures::CouplingCell<dim>*>& inputCells,
                          std::vector<coupling::datastructures::CouplingCell<dim>*>& m2mDomainCells,
                          std::vector<tarch::la::Vector<dim, unsigned int>>& m2mIndices, std::vector<coupling::datastructures::CouplingCell<dim>*>& outerCells,
                          std::vector<tarch::la::Vector<dim, unsigned int>>& outerIndices);

  /** For all d < dim, lowerBoundaries[d] < upperBoundaries[d].
   *  If M2M-domain is not yet defined, this does nothing.
   *
   *  @brief Writes two limiting boundary vectors to the arguments passed.
   *  @param lowerBoundaries Index to write lower boundaries to.
   *  @param upperBoundaries Index to write upper boundaries to.
   * */
  void getMD2MacroDomainBoundariesAllRanks(tarch::la::Vector<dim, unsigned int>& lowerBoundaries, tarch::la::Vector<dim, unsigned int>& upperBoundaries) const {
    if (_lowerBoundaryAllRanks != _upperBoundaryAllRanks) {
      lowerBoundaries = *_lowerBoundaryAllRanks;
      upperBoundaries = *_upperBoundaryAllRanks;
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    else
      std::cout << "WARNING: ICM2M (" << _myRank
                << "): getGlobalMD2MacroDomainBoundaries while domain "
                   "boundaries are unitialized!"
                << std::endl;
#else
    else
      std::cout << "WARNING: ICM2M: getGlobalMD2MacroDomainBoundaries while "
                   "domain boundaries are unitialized!"
                << std::endl;
#endif
  }

  /** For all d < dim, lowerBoundaries[d] < upperBoundaries[d].
   *  If M2M-domain is not yet defined, this does nothing.
   *
   *  @brief Writes two limiting boundary vectors to the arguments passed.
   *  @param lowerBoundaries Index to write lower boundaries to.
   *  @param upperBoundaries Index to write upper boundaries to.
   * */
  void getMD2MacroDomainBoundariesThisRank(tarch::la::Vector<dim, unsigned int>& lowerBoundaries, tarch::la::Vector<dim, unsigned int>& upperBoundaries) const {
    if (_lowerBoundaryThisRank != _upperBoundaryThisRank) {
      lowerBoundaries = *_lowerBoundaryThisRank;
      upperBoundaries = *_upperBoundaryThisRank;
    }
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    else
      std::cout << "WARNING: ICM2M (" << _myRank
                << "): getLocalMD2MacroDomainBoundaries while domain "
                   "boundaries are unitialized!"
                << std::endl;
#else
    else
      std::cout << "WARNING: ICM2M: getLocalMD2MacroDomainBoundaries while "
                   "domain boundaries are unitialized!"
                << std::endl;
#endif
  }

  /** Determines the size of the global MD-To-Macro domain, i.e. the number of
   * cells contained in the domain across all ranks.
   *
   *  @return size of the global MD-To-Macro domain		 */
  tarch::la::Vector<dim, unsigned int> getGlobalMD2MacroDomainSize() const {
    // Since both lower and upper boundaries are inclusive, we need to add one.
    // operator+(int) in tarch::la:Vector would be great here...
    auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
    return *_upperBoundaryAllRanks - *_lowerBoundaryAllRanks + plus_one;
  }

  /** Determines the size of this rank's local MD-To-Macro domain, i.e. the
   * number of cells contained in the domain on this rank only.
   *
   *  @return size of the local MD-To-Macro domain 		 */
  tarch::la::Vector<dim, unsigned int> getLocalMD2MacroDomainSize() const {
    auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
    return *_upperBoundaryThisRank - *_lowerBoundaryThisRank + plus_one;
  }

  /** Determines the global offset of the MD-To-Macro domain relative to the
   * global MD domain.
   *
   *  @return global MD-To-Macro domain offset 		 */
  tarch::la::Vector<dim, double> getGlobalMD2MacroDomainOffset() const {
    tarch::la::Vector<dim, double> offset = _ic->getGlobalMDDomainOffset(); // standard MD offset

    if (_lowerBoundaryAllRanks == nullptr)
      throw std::runtime_error("ERROR: Calling while _lowerBoundaryAllRanks is uninitialized!");

    // offset of md2macro domain relative to MD domain
    for (unsigned int d = 0; d < dim; d++) {
      offset[d] += _ic->getCouplingCellSize()[d] * (*_lowerBoundaryAllRanks)[d]; // offset of md2macro domain
                                                                                    // relative to MD domain
    }

    // std::cout << offset << std::endl << std::endl;
    return offset;
  }

  /** Same as getGlobalVectorCellIndex of base class, but converts in terms of
   * to MD-To-Macro indexing. Has the option to "ignore" Ghost Layer cells, i.e.
   * return INT_MAX for all indices in Ghost Layer. E.g.: Ghost layer at x = 0.
   * Then requesting vector index (x,y,z,..) will return (INT_MAX, y, z,...).
   *
   *  @brief Conversion function: Linear global index -> global vector index
   *  @param globalCellIndex Linear global index to be converted to vector
   *  @param noGL True, iff Ghost Layer indices shall be set to INT_MAX
   *  @return Converted vector index 		 */
  tarch::la::Vector<dim, unsigned int> getGlobalVectorCellIndex(unsigned int globalCellIndex, bool noGL = true) const;

  /** Same as getLocalVectorCellIndex of base class, but converts in terms of to
   * MD-To-Macro indexing. Has the option to "ignore" Ghost Layer cells, i.e.
   * return INT_MAX for all indices in Ghost Layer. E.g.: Ghost layer at x = 0.
   * Then requesting vector index (x,y,z,..) will return (INT_MAX, y, z,...).
   *
   *  @brief Conversion function: Linear local index -> global vector index
   *  @param localCellIndex Linear local index to be converted to vector
   *  @param noGL True, iff Ghost Layer indices shall be set to INT_MAX
   *  @return Converted vector index 		 */
  tarch::la::Vector<dim, unsigned int> getLocalVectorCellIndex(unsigned int localCellIndex, bool noGL = true) const;

  /** In a lot of cases, you want instances of this to be indentical to a
   * coupling::IndexConversion. This way you can access this class' underlying
   * IC element with ease.
   *
   *  @return the IndexConversion this wrapper class is used for */
  const coupling::IndexConversion<dim>* getBaseIC() const { return _ic; }

  /** Same as getLocalCellIndex, but converts in terms of MD-To-Macro indexing.
   *
   *  @brief Conversion function: Local vector index -> Linear local index
   *  @param Local vector index to be converted to linear representation.
   *  @return Converted linear index 		 */
  unsigned int getLocalCellIndex(tarch::la::Vector<dim, unsigned int> localCellIndex) const {
    auto numberCells = getLocalMD2MacroDomainSize();
    unsigned int index = localCellIndex[dim - 1];
    for (int d = dim - 2; d > -1; d--)
      index = numberCells[d] * index + localCellIndex[d];
    return index;
  }

private:
  /** Used within initMD2MacroDomain().
   *  Initilizes _{lower, upper}Boundary{All, This}Ranks (see below). */
  void initGlobalMD2MacroDomainBoundaries();
  /** IndexConversionMD2Macro is a wrapper class for IndexConversion. This is
   * the underlying IndexConversion. */
  const coupling::IndexConversion<dim>* _ic;
  /** Used to determine which cell indices are received by the Macroscopic
   * Solver. */
  coupling::interface::MacroscopicSolverInterface<dim>* _msi;
  /** Indices of the first cell of the MD-To-Macro domain on all ranks.
   *  Initialised during first call of getMD2MacroDomainBoundaries.
   *  Indexing here is relative to the entire MD-Domain (i.e. what is sometimes
   * referred to as MD-Indexing/MaMiCo-Indexing) */
  tarch::la::Vector<dim, unsigned int>* _lowerBoundaryAllRanks;
  /** Indices of the last cell of the MD-To-Macro domain on all ranks.
   *  Initialised during first call of getMD2MacroDomainBoundaries.
   *  Indexing here is relative to the entire MD-Domain (i.e. what is sometimes
   * referred to as MD-Indexing/MaMiCo-Indexing) */
  tarch::la::Vector<dim, unsigned int>* _upperBoundaryAllRanks;
  /** Indices of the first cell of the MD-To-Macro domain on this rank.
   *  Initialised during first call of getMD2MacroDomainBoundaries.
   *  Indexing here is relative to the entire MD-Domain (i.e. what is sometimes
   * referred to as MD-Indexing/MaMiCo-Indexing) */
  tarch::la::Vector<dim, unsigned int>* _lowerBoundaryThisRank;
  /** Indices of the first last cell of the MD-To-Macro domain on this rank.
   *  Initialised during first call of getMD2MacroDomainBoundaries.
   *  Indexing here is relative to the entire MD-Domain (i.e. what is sometimes
   * referred to as MD-Indexing/MaMiCo-Indexing) */
  tarch::la::Vector<dim, unsigned int>* _upperBoundaryThisRank;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  /** mpi communicator */
  const MPI_Comm _comm;
  /** Lowest rank in _comm which contains cell(s).
   *  This rank is assumed to manage cell (0,...,0) of both the MD and the
   * MD-To-Macro-domain. */
  const int _lowestRank;
  /** This ICM2M instance's rank */
  const int _myRank;
#endif
};

#include "IndexConversionMD2Macro.cpph"
