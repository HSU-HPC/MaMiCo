// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_
#define _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_

#include "tarch/la/Vector.h"
#include "tarch/utils/Uncopyable.h"
#include <cstdlib>
#include <vector>
#include "coupling/CouplingMDDefinitions.h"
// parallel topologies
#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

namespace coupling {
  template<unsigned int dim>
  class IndexConversion;
}


/** handles the index convertion from vector to linearised indices as well as from local to global indices.
 *  This class inherits privately from Uncopyable (see CouplingMDDefinitions). We do this in order to hold the
 *  pointer to the ParallelTopology consistent in exactly one instance of IndexConversion (or require to have
 *  valid initialisation via the IndexConversion(...) constructor.
 *
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::IndexConversion: private tarch::utils::Uncopyable {

  public:
    /** constructor for Multi-MD simulations. */
    IndexConversion(
      tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,
      tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int rank,
      tarch::la::Vector<dim,double> globalMDDomainSize,
      tarch::la::Vector<dim,double> globalMDDomainOffset,
      coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
      unsigned int topologyOffset
    );
    /** constructor for single-MD simulations */
    IndexConversion(
      tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,
      tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int rank,
      tarch::la::Vector<dim,double> globalMDDomainSize,
      tarch::la::Vector<dim,double> globalMDDomainOffset,
      coupling::paralleltopology::ParallelTopologyType parallelTopologyType
    ): IndexConversion<dim>(globalNumberMacroscopicCells,numberProcesses,rank,globalMDDomainSize,globalMDDomainOffset,parallelTopologyType,0){}

    ~IndexConversion();


    // ----------------------------- INDEX CONVERSIONS ------------------------------------

    /** returns the linearised global cell index from a global cell index vector. We assume a lexicographic ordering
     *  of the macroscopic cells for this purpose. Besides, we assume the number of global macroscopic cells to be
     *  extended by a global ghost layer (so that we have globalNumberMacroscopicCells_d+2 cells in each direction).
     */
    unsigned int getGlobalCellIndex(tarch::la::Vector<dim,unsigned int> globalCellIndex) const;

    /** returns the linearised global cell index from a local cell index vector. We assume a lexicographic ordering
     *  of the macroscopic cells for this purpose. Besides, we assume the number of local macroscopic cells to be
     *  extended by a global ghost layer (so that we have globalNumberMacroscopicCells_d+2 cells in each direction).
     */
    unsigned int getLocalCellIndex(tarch::la::Vector<dim,unsigned int> localCellIndex) const;

    /** returns the global vector cell index from a linearised index. This function just performs the inverse operation to
     *  getGlobalCellIndex(vectorIndex), i.e. getGlobalVectorCellIndex(getGlobalCellIndex(somevector) ) should again return
     *  somevector.
     */
    tarch::la::Vector<dim,unsigned int> getGlobalVectorCellIndex(unsigned int globalCellIndex) const;

    /** returns the local vector cell index from a linearised index. This function just performs the inverse operation to
     *  getLocalCellIndex(vectorIndex), i.e. getLocalVectorCellIndex(getLocalCellIndex(somevector) ) should again return
     *  somevector.
     */
    tarch::la::Vector<dim,unsigned int> getLocalVectorCellIndex(unsigned int localCellIndex) const;

    /** converts a local vector cell index into a global vector cell index */
    tarch::la::Vector<dim,unsigned int> convertLocalToGlobalVectorCellIndex(
      tarch::la::Vector<dim,unsigned int> localCellIndex
    ) const;

    /** converts a global vector cell index into a local vector cell index */
    tarch::la::Vector<dim,unsigned int> convertGlobalToLocalVectorCellIndex(
      tarch::la::Vector<dim,unsigned int> globalCellIndex
    ) const;

    /** converts a local to global (linearised) cell index. */
    unsigned int convertLocalToGlobalCellIndex(unsigned int localCellIndex) const;

    /** converts a global to a local (linearised) cell index.*/
    unsigned int convertGlobalToLocalCellIndex(unsigned int globalCellIndex) const;


    // ------------------------------ GETTERS ----------------------------------------------

    /** returns the "average" number of local macroscopic cells. "Average" in this case means
     *  the number of cells on all processes except for the processes on the final coordinate range;
     *  we may have more cells there, due to the domain decomposition. The number of ghost cells
     *  is not included in this vector.
     */
    tarch::la::Vector<dim,unsigned int> getAverageLocalNumberMacroscopicCells() const;
    /** returns the local number of macroscopic cells on this process. The ghost layer is not
     *  included in this vector.
     */
    tarch::la::Vector<dim,unsigned int> getLocalNumberMacroscopicCells() const;
    /** returns the global number of macroscopic cells, exluding ghost cells. */
    tarch::la::Vector<dim,unsigned int> getGlobalNumberMacroscopicCells() const;
    /** returns the coordinates of this process. */
    tarch::la::Vector<dim,unsigned int> getThisProcess() const;
    /** returns the rank which corresponds to the linearised index of the vector returned by
     *  getThisProcess().
     */
    unsigned int getThisRank() const;
    /** returns the vector with the number of processes used in all directions. */
    tarch::la::Vector<dim,unsigned int> getNumberProcesses() const;
    /** returns the global domain size of the MD domain (excl. ghost layer which naturally
     *  is not part of MD).
     */
    tarch::la::Vector<dim,double> getGlobalMDDomainSize() const ;
    /** returns the offset, i.e. the lower,left... corner coordinate, of the MD domain.*/
    tarch::la::Vector<dim,double> getGlobalMDDomainOffset() const;
    /** returns the vector size of each macroscopic cell. */
    tarch::la::Vector<dim,double> getMacroscopicCellSize() const;

    coupling::paralleltopology::ParallelTopologyType getParallelTopologyType() const;


    // ---------------------------- GEOMETRY TO CELL INDEX ----------------------------------

    /** returns the global vector cell index of the macroscopic cell which contains the vector "position". */
    tarch::la::Vector<dim,unsigned int> getGlobalVectorCellIndex(tarch::la::Vector<dim,double> position) const;


    // ---------------------------- RANKS AND CELL INDICES ----------------------------------

    /** returns the process coordinates for the rank 'rank'. Forwards the call to the ParallelTopology implementation. */
    tarch::la::Vector<dim,unsigned int> getProcessCoordinates(unsigned int rank) const;
    /** returns the linearised index (i.e. rank) of the given process coordinates. Forwards the call to the
     *  ParallelTopology implementation.
     */
    unsigned int getRank(tarch::la::Vector<dim,unsigned int> processCoordinates) const;

    /** returns the unique rank for a macroscopic cell. If the cell is contained in the MD volume, the rank is uniquely
     *  chosen by using the process which contains the macroscopic cell as non-ghost (i.e. real inner) cell. If the cell
     *  is a global ghost cell, we choose the rank according to the block-decomposition of the grid applied to the
     *  complete grid incl. the ghost layer.
     */
    unsigned int getUniqueRankForMacroscopicCell(tarch::la::Vector<dim,unsigned int> globalCellIndex) const;

    /** returns all ranks which contain a copy (as ghost or non-ghost cell) of the macroscopic cell with this global
     *  cell index. This function only supports 1D, 2D, 3D.
     */
    std::vector<unsigned int> getRanksForMacroscopicCell(tarch::la::Vector<dim,unsigned int> globalCellIndex) const;

    // ------------------------ VALIDATION --------------------------------------------

    /** checks if the global vector cell index is in a well-defined range. Should only be used for debugging
     *  and validation purposes.
     */
    bool isValidGlobalVectorCellIndex(tarch::la::Vector<dim,unsigned int> globalCellIndex) const;
    /** checks if the local vector cell index is in a well-defined range. Should only be used for debugging
     *  and validation purposes.
     */
    bool isValidLocalVectorCellIndex(tarch::la::Vector<dim,unsigned int> localCellIndex) const;
    /** checks if the global cell index is in a well-defined range. Should only be used for debugging
     *  and validation purposes.
     */
    bool isValidGlobalCellIndex(unsigned int globalIndex) const;
    /** checks if the local cell index is in a well-defined range. Should only be used for debugging
     *  and validation purposes.
     */
    bool isValidLocalCellIndex(unsigned int localIndex) const;

  private:
    /** converts a vector cell index into a linearised cell index using a predefined number of cells.
     *  We add the ghost layer influence inside this function.
     *  This function is supposed to be only called with a member variable for numberCells.
     */
    unsigned int getCellIndex(
      tarch::la::Vector<dim,unsigned int> vectorCellIndex, const tarch::la::Vector<dim,unsigned int>& numberCells
    ) const;

    /** returns the average number of local macroscopic cells on each process.
     *  Example: we have 31x46x50 macroscopic cells split onto 10x12x5 processes.
     *  Then, this method will return the floor (cells/processes), which in our case
     *  yields (3x3x10) (-> floor(31/10)=3, floor(46/12)=3, floor(50/5) = 10).
     */
    tarch::la::Vector<dim,unsigned int> initAverageLocalNumberMacroscopicCells(
      tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,
      tarch::la::Vector<dim,unsigned int> numberProcesses
    ) const;

    /** this method returns the exact number of local macroscopic cells.
     *  For Nx x Ny x Nz processes, the method returns the average local number of
     *  macroscopic cells for all processes from [1;Nx-1]x[1;Ny-1]x[1;Nz-1].
     *  For the other processes (at the right/back/top boundary), the number of cells
     *  is filled up so that the global number of macroscopic cells is reached.
     */
    tarch::la::Vector<dim,unsigned int> initLocalNumberMacroscopicCells(
      tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells,
      tarch::la::Vector<dim,unsigned int> numberProcesses,
      unsigned int rank
    ) const;

    /** checks if the index is in a valid range. We add the ghost layer influence inside this function. */
    bool isValidVectorCellIndex(
      tarch::la::Vector<dim,unsigned int> index, const tarch::la::Vector<dim,unsigned int> &range
    ) const;

    /** returns true, if the index is smaller than the range. */
    bool isValidCellIndex(unsigned int index, const unsigned int &range) const;


    /** initialises _macroscopicCellSize using _globalMDDomainSize and _globalNumberMacroscopicCells. */
    tarch::la::Vector<dim,double> initMacroscopicCellSize(
      const tarch::la::Vector<dim,double>& globalMDDomainSize, const tarch::la::Vector<dim,unsigned int>& globalNumberMacroscopicCells
    ) const;


    /** parallel topology used in the simulation */
    const coupling::paralleltopology::ParallelTopologyType _parallelTopologyType;
    const coupling::paralleltopology::ParallelTopology<dim> *_parallelTopology;
    /** the total number of processes. We assume a block-domain decomposition. */
    const tarch::la::Vector<dim,unsigned int> _numberProcesses;
    /** rank of current process */
    const unsigned int _rank;
    /** the coordinates of the current process. */
    const tarch::la::Vector<dim,unsigned int> _thisProcess;

    /** global number of (inner) macroscopic cells. In total this number is extended by 2 in order to account for an
     *  additional ghost layer surrounding the global domain.
     */
    const tarch::la::Vector<dim,unsigned int> _globalNumberMacroscopicCells;
    /** average number of local macroscopic cells */
    const tarch::la::Vector<dim,unsigned int> _averageLocalNumberMacroscopicCells;
    /** local number of macroscopic cells for this process. In total, this number is extended by 2 in order to account
     *  for an additional ghost layer surrounding the local domain.
     */
    const tarch::la::Vector<dim,unsigned int> _localNumberMacroscopicCells;

    /** division factors (in 3D, its (1,g(0)+2,(g(0)+2)*(g(1)+2)) ) for the global number of macroscopic cells. */
    const tarch::la::Vector<dim,unsigned int> _divisionFactor4GlobalNumberMacroscopicCells;

    /** division factors (in 3D, its (1,g(0)+2,(g(0)+2)*(g(1)+2)) ) for the local number of macroscopic cells. */
    const tarch::la::Vector<dim,unsigned int> _divisionFactor4LocalNumberMacroscopicCells;

    /** global MD domain size. */
    const tarch::la::Vector<dim,double> _globalMDDomainSize;
    /** offset (lower,left... corner) of MD domain. */
    const tarch::la::Vector<dim,double> _globalMDDomainOffset;
    /** size of macroscopic cell. */
    const tarch::la::Vector<dim,double> _macroscopicCellSize;
};

#include "IndexConversion.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_INDEXCONVERSION_H_
