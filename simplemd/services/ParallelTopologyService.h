// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_SERVICES_PARALLELTOPOLOGYSERVICE_H_
#define _MOLECULARDYNAMICS_SERVICES_PARALLELTOPOLOGYSERVICE_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/la/Vector.h"
#if (MD_PARALLEL == MD_YES)
#include <mpi.h>
#endif
#include "simplemd/LinkedCell.h"
#include "simplemd/Molecule.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/services/ParallelAndLocalBufferService.h"
#include <cstdlib>
#include <iostream>
#include <vector>

namespace simplemd {
namespace services {
class ParallelTopologyService;

// forward declarations to circumvent circular includes
class MoleculeService;
class LinkedCellService;
} // namespace services
} // namespace simplemd

/** services managing the distributed memory parallelisation things.
 *  It is based on a simple domain decomposition method for the molecules.
 *  Therefore, if we have N=Nx x Ny x Nz processes available, we enumerate
 *  the processes lexicographically in their ranks and so set up an easy-to-use
 *  topology of our parallel simulation.
 *
 *  @author Philipp Neumann, Nikola Tchipev
 */
class simplemd::services::ParallelTopologyService {
public:
  /** initialise the service with global domain size and offset,
   *  predefined meshwidth for the linked cells (which is to be adopted such
   * that
   *  it fits to the local process properties), the number of processes to be
   * used
   *  in each spatial direction, the rank of the current process, the global
   *  description of the outer boundaries, the number of molecules per direction
   *  and initialises the ParallelAndLocalBufferService.
   *
   *  Needs to be called before the local services are initialised.
   */
  ParallelTopologyService(const tarch::la::Vector<MD_DIM, double> &domainSize, const tarch::la::Vector<MD_DIM, double> &domainOffset,
                          const tarch::la::Vector<MD_DIM, double> &meshWidth, const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
                          const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &boundary
#if (MD_PARALLEL == MD_YES)
                          ,
                          MPI_Comm communicator = MPI_COMM_WORLD
#endif
  );
  ~ParallelTopologyService() {}

  /** finish initialization by allocating buffers
   *
   *  Needs to be called after MoleculeService has been initialised.
   */
  void initBuffers(const unsigned int &localNumberOfMolecules);

  /** shuts down the buffer service and then the parallel service */
  void shutdown();

  /** returns the rank of the current process */
  int getRank() const;

  /** returns the process coordinates of the local process */
  const tarch::la::Vector<MD_DIM, unsigned int> &getProcessCoordinates() const;

  /** returns local information on boundary relations */
  const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &getLocalBoundaryInformation() const { return _boundary; }

  /** returns the mesh width for the linked cells in the current simulation.
   * Depending on the meshwidth
   *  handed over during initialisation, this meshwidth might be different; we
   * try to choose it as close
   *  to the meshwidth defined by the user, but big enough so that we have the
   * same (integer) number of
   *  cells on each process. So, _meshWidth typicall is bigger or equal the
   * user-defined meshwidth.
   *  This may result in higher computation times! However, the physics remain
   * the same as the cutoff-radius
   *  is not affected by these changes.
   */
  const tarch::la::Vector<MD_DIM, double> &getMeshWidth() const { return _meshWidth; }

  /** returns the global number of cells in each spatial direction */
  const tarch::la::Vector<MD_DIM, unsigned int> &getGlobalNumberOfCells() const { return _globalNumberOfCells; }

  /** returns the local number of cells in each spatial direction (i.e. only the
   * cells of this process) */
  const tarch::la::Vector<MD_DIM, unsigned int> &getLocalNumberOfCells() const { return _localNumberOfCells; }

  /** returns the number of processes used in each spatial direction */
  const tarch::la::Vector<MD_DIM, unsigned int> &getNumberOfProcesses() const { return _numberProcesses; }

  /** returns the global domain size*/
  const tarch::la::Vector<MD_DIM, double> &getGlobalDomainSize() const { return _domainSize; }

  /** returns the global domain offset */
  const tarch::la::Vector<MD_DIM, double> &getGlobalDomainOffset() const { return _domainOffset; }

  /** returns the global index of the first (non-ghost) cell */
  const tarch::la::Vector<MD_DIM, unsigned int> &getGlobalIndexOfFirstCell() const { return _globalIndexOfFirstCell; }

  /** returns true, if this process does not carry any work. This can be the
   * case, if we have more ranks available
   *  in the NodePool than specified in the xml config.
   */
  bool isIdle() const;

  /** broadcasts all molecules contained in the (inner) cell "cell" with index
   * cellIndex to all
   * ghost cells that can be identified with this inner cell. Returns a vector
   * containing the local
   * cell coordinates of all neighbored ghost cells that are a local,
   * non-parallel boundary.
   * This is important for the case that we have 1 processor in one direction
   * only. Then, we might
   * have a periodic boundary that needs to be handled locally!
   *
   * This method is also used for handling process-leaving particles in the
   * optimised case when
   * they are sent together with molecules for ghost layers.
   * The local buffer of the buffer service is filled only with process-leaving
   * particles that need
   * to be handled locally.
   * Ghost particles that need to be handled locally are handled via the
   * indices in the returned vector.
   *
   * In the parallel case, this method makes use of the communication buffers:
   * molecules are not
   * messaged to neighbour one by one, but rather pushed into buffers, which
   * store all molecules to be sent to a process.
   *
   * The MPI send calls are executed later (from within BoundaryTreatment).
   */
  std::vector<tarch::la::Vector<MD_DIM, unsigned int>> broadcastInnerCellViaBuffer(LinkedCell &cell, const unsigned int &cellIndex,
                                                                                   const simplemd::services::LinkedCellService &linkedCellService);

  /** sends all molecules from cell cellIndex to the respective neighbouring
   * process. The cell
   * cellIndex needs to be a ghost cell.
   * The function returns true, if the respective cell where the molecules need
   * to be sorted in
   * is part of another process; in this case the molecules are pushed to the
   * send buffers and erased.
   * Otherwise, the other (non-ghost) cell is part of the same process. In this
   * case, no send-operation is triggered and false
   * is returned.
   *
   * The MPI send calls are executed later (from within BoundaryTreatment).
   */
  bool reduceGhostCellViaBuffer(LinkedCell &cell, const unsigned int &cellIndex, const simplemd::services::LinkedCellService &linkedCellService);

  /** unpack and resort local buffer
   */
  void unpackLocalBuffer(simplemd::services::MoleculeService &moleculeService, simplemd::services::LinkedCellService &linkedCellService);

  /** Communication schedule:
   * 1. Irecv on buffers
   * 2. Isend on buffers
   * 3. Wait for all requests to be completed
   * 4. Unpack receive buffers and sort into linked cells
   *
   * For large messages, waiting for sender requests is also a must-have.
   */
  void communicationSteps_1_2();

  /** See comment of communicationSteps_1_2() */
  void communicationSteps_3_4(simplemd::services::MoleculeService &moleculeService, simplemd::services::LinkedCellService &linkedCellService);

  /** Compute (non-overlapping) intersection of a global region of interest
      (ROI) with local domain.
      For example for purposes of profile plotter. */
  bool globalToLocalRegionOfInterest(const tarch::la::Vector<MD_DIM, unsigned int> &globalStartCell, const tarch::la::Vector<MD_DIM, unsigned int> &globalRange,
                                     tarch::la::Vector<MD_DIM, unsigned int> &localStartCell, tarch::la::Vector<MD_DIM, unsigned int> &localRange) const;

  /** Compute global cell index from local cell index in vector form */
  tarch::la::Vector<MD_DIM, unsigned int> localToGlobalCellIndexVector(const tarch::la::Vector<MD_DIM, unsigned int> &localCellIndexVector) const;

#if (MD_PARALLEL == MD_YES)
  int MD_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op) {
    return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, _communicator);
  }
#endif

private:
  /** computes all neighbour ranks and stores the results in neighbourRanks.
   * Computes _neighbourRanksUnique, _numUniqueNeighbours and cells per buffer
   * -
   * values needed by the buffer service.
   */
  void createNeighbourRanks(const tarch::la::Vector<MD_DIM, unsigned int> &processCoordinates,
                            const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &localBoundary,
                            const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses, std::vector<int> &neighbourRanks,
                            std::vector<int> &neighbourRanksUnique, unsigned int &numUniqueNeighbours, unsigned int numberOfCellsPerBuffer[]) const;

  /** if addedNeighbour is not already present in neighbourRanksUnique,
   * add it, also incrementing numUniqueNeighbours.
   */
  void addNeighbourToNeighbourRanksUnique(std::vector<int> &neighbourRanksUnique, unsigned int &numUniqueNeighbours, const int &addedNeighbour) const;

  /** computes the actual mesh size that is used in the simulation. We make the
   * width as similar to the prescribed mesh width "prescribedWidth" as possible
   * so to
   *  have an integer number of grid cells.
   */
  tarch::la::Vector<MD_DIM, double> computeMeshwidth(const tarch::la::Vector<MD_DIM, double> &prescribedWidth,
                                                     const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
                                                     const tarch::la::Vector<MD_DIM, double> &domainSize) const;

  /** computes the number of cells in each process block. The
   *  number of cells computed by this service should also be used
   *  by the LinkedCellService!
   */
  tarch::la::Vector<MD_DIM, unsigned int> computeNumberOfCells(const tarch::la::Vector<MD_DIM, double> &meshWidth,
                                                               const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
                                                               const tarch::la::Vector<MD_DIM, double> &domainSize) const;

  /** computes the global number of grid cells */
  tarch::la::Vector<MD_DIM, unsigned int> computeGlobalNumberOfCells(const tarch::la::Vector<MD_DIM, double> &meshWidth,
                                                                     const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
                                                                     const tarch::la::Vector<MD_DIM, double> &domainSize) const;

  /** returns true if neighbourRank is present in the _neighbourRanks vector.
   * Returns false otherwise. */
  bool isParallelNeighbour(const int &neighbourRank) const;

  /** returns the local information needed for periodic boundary treatment. This
   * function
   *  is called within the init() call and initialises the _periodicBoundary
   * field.
   */
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>
  computeLocalBoundaryInformation(const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &boundary,
                                  const tarch::la::Vector<MD_DIM, unsigned int> &processCoordinates,
                                  const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses);

  /** returns the number of transferred cells, if the neighbour is placed in
   * distance x,y,z
   *  from the current cell.
   */
  unsigned int getNumberOfTransferredCells(const int &x
#if (MD_DIM > 1)
                                           ,
                                           const int &y
#endif
#if (MD_DIM > 2)
                                           ,
                                           const int &z
#endif
  ) const;

  /** corrects the position vector 'position' according to periodic boundary
   * conditions.
   *  If the boundaryType is PERIODIC_BOUNDARY the index neighbourIndex is used
   * to determine
   *  the location of the present boundary and to modify the position vector
   * accordingly.
   */
  void adaptPositionForPeriodicBoundaries(tarch::la::Vector<MD_DIM, double> &position, const simplemd::BoundaryType &boundaryType, const int &x
#if (MD_DIM > 1)
                                          ,
                                          const int &y
#endif
#if (MD_DIM > 2)
                                          ,
                                          const int &z
#endif
  ) const;

  /** get buffer index corresponding to neighbourRank from _neighbourRanksUnique
   */
  unsigned int getCurrentBufferIndexFromNeighbourRank(const int &neighbourRank) const;

  /** Read off all molecules from buffer and resort them in the respective
   * linked cells.
   */
  void unpackBuffer(ParallelAndLocalBufferService::SimpleBuffer *buf, simplemd::services::MoleculeService &moleculeService,
                    simplemd::services::LinkedCellService &linkedCellService);

  /** place position, velocity, forceOld and isFixed at the end of the
   * respective local buffer.
   */
  void pushMoleculeToLocalBuffer(const Molecule *mol, const tarch::la::Vector<MD_DIM, double> &pos);

#if (MD_PARALLEL == MD_YES)
  /** place position, velocity, forceOld and isFixed at the end of the
   * respective sendBuffer.
   */
  void pushMoleculeToSendBuffer(const unsigned int &bufferIndex, const Molecule *mol, const tarch::la::Vector<MD_DIM, double> &pos);

  /** Issues an Isend call on buffer.
   * Does not wait for the request to be fulfilled, but uses up the respective
   * request
   */
  void bufferIsend(ParallelAndLocalBufferService::SimpleBuffer *buffer, const int &neighbourRank, MPI_Request &request) const;

  /** Issues an Irecv call on buffer.
   * Does not wait for the request to be fulfilled, but uses up the respective
   * request
   */
  void bufferIrecv(ParallelAndLocalBufferService::SimpleBuffer *buffer, const int &neighbourRank, MPI_Request &request) const;
#endif

  /** domain size */
  const tarch::la::Vector<MD_DIM, double> _domainSize;

  /** domain offset */
  const tarch::la::Vector<MD_DIM, double> _domainOffset;

  const tarch::la::Vector<MD_DIM, double> _meshWidth;

  /** index of first cell w.r.t. to global grid
   */
  tarch::la::Vector<MD_DIM, unsigned int> _globalIndexOfFirstCell;

  /** definition of the process matrix, i.e. how many processes work in each
   *  spatial direction.
   */
  const tarch::la::Vector<MD_DIM, unsigned int> _numberProcesses;

  /** local number of cells (i.e. on each process) in all d directions */
  const tarch::la::Vector<MD_DIM, unsigned int> _localNumberOfCells;

  /** global number of cells in all d directions */
  const tarch::la::Vector<MD_DIM, unsigned int> _globalNumberOfCells;

  /** rank of this process */
  int _rank;

  /** boundary information on the local process */
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> _boundary;

  /** coordinates of this process in the process matrix. This vector is the same
   * as
   *  _rank, but written in vector form.
   */
  tarch::la::Vector<MD_DIM, unsigned int> _processCoordinates;

  /** vector containing all ranks of the neighbouring processes. This includes
   *  periodic neighbour relations.
   */
  std::vector<int> _neighbourRanks;

  /** Buffer service exported in a separate class
   * @see ParallelAndLocalBufferService
   */
  ParallelAndLocalBufferService _bufferService;

  /** like _neighbourRanks, but without possible repetitions
   * @see _neighbourRanks
   */
  std::vector<int> _neighbourRanksUnique;

#if (MD_PARALLEL == MD_YES)
  /** tag for communicating buffers */
  int _bufferTag;
#endif

  /** number of processes current process will communicate with
   * also equal to number of used buffers and MPI requests
   */
  unsigned int _numUniqueNeighbours;

  /** how many cells a buffer is responsible for */
  unsigned int _numberOfCellsPerBuffer[MD_LINKED_CELL_NEIGHBOURS];

#if (MD_PARALLEL == MD_YES)
  MPI_Comm _communicator;
  /** requests for sending and receiving buffers */
  MPI_Request _receiveRequests[MD_LINKED_CELL_NEIGHBOURS];
  MPI_Request _sendRequests[MD_LINKED_CELL_NEIGHBOURS];
#endif
};
#endif // _MOLECULARDYNAMICS_SERVICES_PARALLELTOPOLOGYSERVICE_H_
