// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_
#define _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_

#include "coupling/CouplingMDDefinitions.h"
#include "tarch/la/Vector.h"
#include "tarch/utils/Uncopyable.h"
#include <climits>
#include <cstdlib>
#include <vector>
// parallel topologies
#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

namespace coupling {
template <unsigned int dim> class IndexConversion;
}

/** Handles the index conversion from vector to linearised indices as well as from local to global
 * indices. This class inherits privately from Uncopyable (see
 * CouplingMDDefinitions). We do this in order to hold the pointer to the
 * ParallelTopology consistent in exactly one instance of IndexConversion (or
 * require to have valid initialisation via the IndexConversion(...)
 * constructor.
 *  @brief handles the index conversion
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::IndexConversion : private tarch::utils::Uncopyable {
public:
  /** @brief constructor for Multi-MD simulations.
   *  @param globalNumberCouplingCells total number of coupling cells of
   * the md simulation
   *  @param numberProcesses number of mpi processes in application
   *  @param rank rank of this mpi process
   *  @param globalMDDomainSize size of the total md domain per direction
   *  @param parallelTopologyType the type of parallel topology see
   * ParallelTopologyFactory.h
   *  @param topologyOffset the total number of processes for a md simulation,
   * e.g. for two md simulations on eight ranks it would be 4; */
  IndexConversion(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells, tarch::la::Vector<dim, unsigned int> numberProcesses, unsigned int rank,
                  tarch::la::Vector<dim, double> globalMDDomainSize, tarch::la::Vector<dim, double> globalMDDomainOffset,
                  coupling::paralleltopology::ParallelTopologyType parallelTopologyType, unsigned int topologyOffset);
  /** @brief constructor for single-MD simulations
   *  @param globalNumberCouplingCells total number of coupling cells of
   * the md simulation
   *  @param numberProcesses number of mpi processes in application
   *  @param rank rank of this mpi process
   *  @param globalMDDomainSize size of the total md domain per diretion
   *  @param globalMDDomainOffset offset of the md domain compared to the
   * coordinate origin
   *  @param parallelTopologyType the type of parallel topology see
   * ParallelTopologyFactory.h */
  IndexConversion(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells, tarch::la::Vector<dim, unsigned int> numberProcesses, unsigned int rank,
                  tarch::la::Vector<dim, double> globalMDDomainSize, tarch::la::Vector<dim, double> globalMDDomainOffset,
                  coupling::paralleltopology::ParallelTopologyType parallelTopologyType)
      : IndexConversion<dim>(globalNumberCouplingCells, numberProcesses, rank, globalMDDomainSize, globalMDDomainOffset, parallelTopologyType, 0) {}
  /** @brief Destructor */
  ~IndexConversion();

  // ----------------------------- INDEX CONVERSIONS
  // ------------------------------------
  /** We assume a lexicographic ordering of the coupling cells for this
   * purpose. Besides, we assume the number of global coupling cells to be
   * extended by a global ghost layer (so that we have
   * globalNumberCouplingCells_d+2 cells in each direction).
   *  @brief returns the linearised global cell index from a global cell index
   * vector.
   *  @param globalCellIndex the global cell index as a vector (according to the
   * dimension of the simulation)
   *  @returns the linearised global cell index for the given dimensional index
   */
  unsigned int getGlobalCellIndex(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;

  /** We assume a lexicographic ordering of the coupling cells for this
   * purpose. Besides, we assume the number of local coupling cells to be
   * extended by a global ghost layer (so that we have
   * globalNumberCouplingCells_d+2 cells in each direction).
   *  @brief returns the linearised global cell index from a local cell index
   * vector.
   *  @param localCellIndex the local (only valid on this rank) cell index as a
   * vector
   *  @returns the linearised local cell index for the given dimensional index
   */
  unsigned int getLocalCellIndex(tarch::la::Vector<dim, unsigned int> localCellIndex) const;

  /** This function just performs the inverse operation to
   * getGlobalCellIndex(vectorIndex), i.e.
   * getGlobalVectorCellIndex(getGlobalCellIndex(somevector) ) should again
   * return somevector.
   *  @brief returns the global vector cell index from a linearised index.
   *  @param globalIndex global continuous index of a cell
   *  @returns the dimensional global index for the given linear index */
  tarch::la::Vector<dim, unsigned int> getGlobalVectorCellIndex(unsigned int globalCellIndex) const;

  /** Does the same as getGlobalVectorCellIndex() but
   *    * non-inner cells (i.e. outside of the MD domain) are converted to
   * MAX_INT
   *    * (0, ... , 0) is the first cell in the MD domain etc.
   */
  tarch::la::Vector<dim, unsigned int> getGlobalInnerVectorCellIndex(unsigned int globalCellIndex) const;

  // TODO: unused function
  /** Translates index vector of inner cell to make it usable in global
   * context*/
  tarch::la::Vector<dim, unsigned int> convertInnerVectorCellIndexToGlobal(tarch::la::Vector<dim, unsigned int> cellIndexVectorMDContext) const;

  /** This function just performs the inverse operation to
   * getLocalCellIndex(vectorIndex), i.e.
   * getLocalVectorCellIndex(getLocalCellIndex(somevector) ) should again return
   * somevector.
   *  @brief returns the local vector cell index from a linearised index.
   *  @param localCellIndex local (only valid on this rank) continuous index of
   * a cell
   *  @returns the dimensional global index for the given linear index */
  tarch::la::Vector<dim, unsigned int> getLocalVectorCellIndex(unsigned int localCellIndex) const;

  /** @brief converts a local vector cell index into a global vector cell index
   *  @param localCellIndex the local (only valid on this rank) cell index as a
   * vector
   *  @returns the global index for the given local index  */
  tarch::la::Vector<dim, unsigned int> convertLocalToGlobalVectorCellIndex(tarch::la::Vector<dim, unsigned int> localCellIndex) const;

  /** @brief converts a global vector cell index into a local vector cell index
   *  @param globalCellIndex the global cell index as a vector (according to the
   * dimension of the simulation)
   *  @returns the local index for the given global index */
  tarch::la::Vector<dim, unsigned int> convertGlobalToLocalVectorCellIndex(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;

  /** @brief converts a local to global (linearised) cell index.
   *  @param localCellIndex local (only valid on this rank) continuous index of
   * a cell
   *  @returns the local index for the given global index */
  unsigned int convertLocalToGlobalCellIndex(unsigned int localCellIndex) const;

  /** @brief converts a global to a local (linearised) cell index.
   *  @param globalIndex global continuous index of a cell
   *  @returns the global index for the given local index  */
  unsigned int convertGlobalToLocalCellIndex(unsigned int globalCellIndex) const;

  // ------------------------------ GETTERS
  // ----------------------------------------------
  /** "Average" in this case means the number of cells on all processes except
   * for the processes on the final coordinate range; we may have more cells
   * there, due to the domain decomposition. The number of ghost cells is not
   * included in this vector.
   *  @brief returns the "average" number of local coupling cells.
   *  @returns the average number of coupling cells per rank (dimensional)*/
  tarch::la::Vector<dim, unsigned int> getAverageLocalNumberCouplingCells() const;

  /** The ghost layer is not included in this vector.
   *  @brief returns the local number of coupling cells on this process.
   *  @returns the number of coupling cells on this rank (dimensional) */
  tarch::la::Vector<dim, unsigned int> getLocalNumberCouplingCells() const;
  /** @brief returns the global number of coupling cells, exluding ghost
   * cells.
   *  @returns the total number of coupling cells over all ranks
   * (dimensional) */
  tarch::la::Vector<dim, unsigned int> getGlobalNumberCouplingCells() const;
  /** @brief returns the coordinates of this process.
   *  @returns the coordinates of this process */
  tarch::la::Vector<dim, unsigned int> getThisProcess() const;
  /** @brief returns the rank which corresponds to the linearised index of the
   * vector returned by getThisProcess().
   *  @returns the linearised index of this process's coordinates */
  unsigned int getThisRank() const;
  /** @brief returns the vector with the number of processes used in all
   * directions.
   *  @returns the total number of processes per spacial direction  */
  tarch::la::Vector<dim, unsigned int> getNumberProcesses() const;
  /** returns the global domain size of the MD domain (excl. ghost layer which
   * naturally is not part of MD).
   *  @returns the total size of the md simulation domain (dimensional) */
  tarch::la::Vector<dim, double> getGlobalMDDomainSize() const;
  /** @brief returns the offset, i.e. the lower,left... corner coordinate, of
   * the MD domain.
   *  @returns the offset of the MD domain */
  tarch::la::Vector<dim, double> getGlobalMDDomainOffset() const;
  /** @brief returns the vector size of each coupling cell.
   *  @returns the size of the coupling cells (dimensional) */
  tarch::la::Vector<dim, double> getCouplingCellSize() const;

  coupling::paralleltopology::ParallelTopologyType getParallelTopologyType() const;

  // ---------------------------- GEOMETRY TO CELL INDEX
  // ----------------------------------
  /** @brief returns the global vector cell index of the coupling cell in
   * which a point is
   *  @param position vector of the point to get the corresponding macroscopic
   * cell
   *  @returns the corresponding index for the global vector cell */
  tarch::la::Vector<dim, unsigned int> getGlobalVectorCellIndex(tarch::la::Vector<dim, double> position) const;

  // ---------------------------- RANKS AND CELL INDICES
  // ----------------------------------
  /** Forwards the call to the ParallelTopology implementation.
   *  @brief returns the process coordinates for the rank 'rank'.
   *  @param rank the linearised/continuous rank of a process
   *  @returns the process coordinates */
  tarch::la::Vector<dim, unsigned int> getProcessCoordinates(unsigned int rank) const;
  /** Forwards the call to the ParallelTopology implementation modulo the number of processes on each dimension (wrap around).
   *  @brief returns the process coordinates for the rank 'rank'.
   *  @param rank the linearised/continuous rank of a process
   *  @returns the process coordinates */
  tarch::la::Vector<dim, unsigned int> getThisProcessCoordinates(unsigned int rank) const;
  /** Forwards the call to the ParallelTopology implementation.
   *  @brief returns the linearised index (i.e. rank) of the given process
   * coordinates.
   *  @param processCoordinates the vector coordinates of a process
   *  @returns the rank for the given coordinates of a process  */
  unsigned int getRank(tarch::la::Vector<dim, unsigned int> processCoordinates) const;
  /** If the cell is contained in the MD volume, the rank is uniquely chosen by
   * using the process which contains the coupling cell as non-ghost (i.e.
   * real inner) cell. If the cell is a global ghost cell, we choose the rank
   * according to the block-decomposition of the grid applied to the complete
   * grid incl. the ghost layer.
   *  @brief returns the unique rank for a coupling cell.
   *  @param globalCellIndex the global vector coordinates of a coupling cell
   *  @returns the rank for the given coupling cell */
  unsigned int getUniqueRankForCouplingCell(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;
  /** returns all ranks which contain a copy (as ghost or non-ghost cell) of the
   * coupling cell with this global cell index. This function only supports
   * 1D, 2D, 3D.
   *  @param globalCellIndex the global vector coordinates of a coupling cell
   *  @returns all ranks for the given coupling cell */
  std::vector<unsigned int> getRanksForCouplingCell(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;

  // ------------------------ VALIDATION
  // --------------------------------------------
  /** Should only be used for debugging and validation purposes.
   *  @brief checks if the global vector cell index is in a well-defined range.
   *  @param globalCellIndex the global vector coordinates of a coupling cell
   *  @returns a bool, indicating if the given global index is valid (true) or
   * not (false) */
  bool isValidGlobalVectorCellIndex(tarch::la::Vector<dim, unsigned int> globalCellIndex) const;
  /** Should only be used for debugging and validation purposes.
   *  @brief checks if the local vector cell index is in a well-defined range.
   *  @param localCellIndex the local (only valid on this rank) vector
   * coordinates of a coupling cell
   *  @returns a bool, indicating if the given local index is valid (true) or
   * not (false)  */
  bool isValidLocalVectorCellIndex(tarch::la::Vector<dim, unsigned int> localCellIndex) const;
  /** Should only be used for debugging and validation purposes.
   *  @brief checks if the global cell index is in a well-defined range.
   *  @param globalIndex global continuous/linearised index of a cell
   *  @returns a bool, indicating if the given global index is valid (true) or
   * not (false) */
  bool isValidGlobalCellIndex(unsigned int globalIndex) const;
  /** Should only be used for debugging and validation purposes.
   *  @brief checks if the local cell index is in a well-defined range.
   *  @param localIndex local continuous/linearised index of a cell
   *  @returns a bool, indicating if the given local index is valid (true) or
   * not (false) */
  bool isValidLocalCellIndex(unsigned int localIndex) const;

private:
  /** We add the ghost layer influence inside this function.
   *  This function is supposed to be only called with a member variable for
   * numberCells.
   *  @brief converts a vector cell index into a linearised cell index using a
   * predefined number of cells.
   *  @param vectorCellIndex vector cell index
   *  @param numberCells total number of cells
   *  @returns a linearised cell index  */
  unsigned int getCellIndex(tarch::la::Vector<dim, unsigned int> vectorCellIndex, const tarch::la::Vector<dim, unsigned int>& numberCells) const;

  /** Example: we have 31x46x50 coupling cells split onto 10x12x5 processes.
   *  Then, this method will return the floor (cells/processes), which in our
   * case yields (3x3x10) (-> floor(31/10)=3, floor(46/12)=3, floor(50/5) = 10).
   *  @brief returns the average number of local coupling cells on each
   * process.
   *  @param initAverageLocalNumberCouplingCells average number of
   * coupling cells per process
   *  @param globalNumberCouplingCells total number of coupling cells
   *  @param numberProcesses total number of mpi processes
   *  @returns the average number of coupling cells on each process */
  tarch::la::Vector<dim, unsigned int> initAverageLocalNumberCouplingCells(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells,
                                                                           tarch::la::Vector<dim, unsigned int> numberProcesses) const;

  /** For Nx x Ny x Nz processes, the method returns the average local number of
   *  coupling cells for all processes from [1;Nx-1]x[1;Ny-1]x[1;Nz-1].
   *  For the other processes (at the right/back/top boundary), the number of
   * cells is filled up so that the global number of coupling cells is
   * reached.
   *  @brief this method returns the exact number of local coupling cells.
   *  @param initLocalNumberCouplingCells the total number of macroscopic
   * cells on this process
   *  @param globalNumberCouplingCells the total number of coupling cells
   *  @param numberProcesses total number of mpi processes
   *  @returns the number of coupling cells on this process */
  tarch::la::Vector<dim, unsigned int> initLocalNumberCouplingCells(tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells,
                                                                    tarch::la::Vector<dim, unsigned int> numberProcesses, unsigned int rank) const;

  /** We add the ghost layer influence inside this function.
   *  @brief checks if the index is in a valid range.
   *  @param index vector index of a cell
   *  @param range vector range of cells
   *  @returns a bool indicating if index is valid (true) or not (false) */
  bool isValidVectorCellIndex(tarch::la::Vector<dim, unsigned int> index, const tarch::la::Vector<dim, unsigned int>& range) const;

  /** @brief returns true, if the index is smaller than the range.
   *  @param index continuous/linearised index of a cell
   *  @param range continuous/linearised range of cells
   *  @returns a bool indicating if index is valid (true) or not (false) */
  bool isValidCellIndex(unsigned int index, const unsigned int& range) const;

  /** @brief initialises the size of coupling cells
   *  @param globalMDDomainSize the total size of the md domain as a vector
   *  @param globalNumberCouplingCells the global number of coupling cells
   *  @returns the coupling cell size for the setup */
  tarch::la::Vector<dim, double> initCouplingCellSize(const tarch::la::Vector<dim, double>& globalMDDomainSize,
                                                      const tarch::la::Vector<dim, unsigned int>& globalNumberCouplingCells) const;

  /** type of parallel topology used in the simulation */
  const coupling::paralleltopology::ParallelTopologyType _parallelTopologyType;
  /** a pointer to the parallel topology */
  const coupling::paralleltopology::ParallelTopology<dim>* _parallelTopology;
  /** the total number of processes. We assume a block-domain decomposition. */
  const tarch::la::Vector<dim, unsigned int> _numberProcesses;
  /** rank of current process */
  const unsigned int _rank;
  /** the coordinates of the current process. */
  const tarch::la::Vector<dim, unsigned int> _thisProcess;
  /** global number of (inner) coupling cells. In total this number is
   * extended by 2 in order to account for an
   *  additional ghost layer surrounding the global domain.  */
  const tarch::la::Vector<dim, unsigned int> _globalNumberCouplingCells;
  /** average number of local coupling cells */
  const tarch::la::Vector<dim, unsigned int> _averageLocalNumberCouplingCells;
  /** local number of coupling cells for this process. In total, this number
   * is extended by 2 in order to account for an additional ghost layer
   * surrounding the local domain. */
  const tarch::la::Vector<dim, unsigned int> _localNumberCouplingCells;
  /** division factors (in 3D, its (1,g(0)+2,(g(0)+2)*(g(1)+2)) ) for the global
   * number of coupling cells. */
  const tarch::la::Vector<dim, unsigned int> _divisionFactor4GlobalNumberCouplingCells;
  /** division factors (in 3D, its (1,g(0)+2,(g(0)+2)*(g(1)+2)) ) for the local
   * number of coupling cells. */
  const tarch::la::Vector<dim, unsigned int> _divisionFactor4LocalNumberCouplingCells;
  /** global MD domain size. */
  const tarch::la::Vector<dim, double> _globalMDDomainSize;
  /** offset (lower,left... corner) of MD domain. */
  const tarch::la::Vector<dim, double> _globalMDDomainOffset;
  /** size of coupling cell. */
  const tarch::la::Vector<dim, double> _couplingCellSize;
};

#include "IndexConversion.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_
