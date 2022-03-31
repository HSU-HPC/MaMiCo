// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/services/ParallelTopologyService.h"

#include "simplemd/services/LinkedCellService.h"
#include "simplemd/services/MoleculeService.h"

simplemd::services::ParallelTopologyService::ParallelTopologyService(
    const tarch::la::Vector<MD_DIM, double> &domainSize,
    const tarch::la::Vector<MD_DIM, double> &domainOffset,
    const tarch::la::Vector<MD_DIM, double> &meshWidth,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &
        boundary
#if (MD_PARALLEL == MD_YES)
    ,
    MPI_Comm communicator
#endif
    )
    : _domainSize(domainSize), _domainOffset(domainOffset),
      _meshWidth(computeMeshwidth(meshWidth, numberProcesses, domainSize)),
      _numberProcesses(numberProcesses),
      _localNumberOfCells(computeNumberOfCells(meshWidth, numberProcesses,
                                               domainSize)),
      _globalNumberOfCells(computeGlobalNumberOfCells(meshWidth,
                                                      numberProcesses,
                                                      domainSize))
#if (MD_PARALLEL == MD_YES)
      ,
      _communicator(communicator)
#endif
      {
  unsigned int help;
// set values
#if (MD_PARALLEL == MD_YES)
  MPI_Comm_rank(_communicator, &_rank);
#else
  _rank = 0;
#endif
  if (_rank < 0) {
    std::cout
        << "ERROR simplemd::services::ParallelTopologyService::init: Rank < 0!"
        << std::endl;
    exit(EXIT_FAILURE);
  }

  // compute process coordinates
  help = (unsigned int) _rank;
#if (MD_DIM > 2)
  _processCoordinates[2] = help / (_numberProcesses[0] * _numberProcesses[1]);
  help = help -
         _processCoordinates[2] * (_numberProcesses[0] * _numberProcesses[1]);
#endif
#if (MD_DIM > 1)
  _processCoordinates[1] = help / _numberProcesses[0];
  help = help - _processCoordinates[1] * _numberProcesses[0];
#endif
  _processCoordinates[0] = help;
  // determine local information on periodic boundary
  _boundary = computeLocalBoundaryInformation(boundary, _processCoordinates,
                                              _numberProcesses);
  // determine ranks of possible communication neighbours
  createNeighbourRanks(_processCoordinates, _boundary, _numberProcesses,
                       _neighbourRanks, _neighbourRanksUnique,
                       _numUniqueNeighbours, _numberOfCellsPerBuffer);
  // compute global index of first cell
  for (unsigned int d = 0; d < MD_DIM; d++) {
    _globalIndexOfFirstCell[d] =
        _processCoordinates[d] * _localNumberOfCells[d];
  }

#if (MD_PARALLEL == MD_YES)
  _bufferTag = 1;
  if (!isIdle()) {
  } else {
    std::cout << "ERROR simplemd::services::ParallelTopologyService::init(): "
                 "Currently no idle processes are allowed!" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

#if (MD_DEBUG == MD_YES)
  std::cout << "ParallelTopologyService: Rank=" << _rank
            << ", process coords: " << _processCoordinates << std::endl;
  std::cout << "No. processes: " << _numberProcesses
            << ", boundary: " << _boundary << std::endl;
  std::cout << "Parallel neighbours:";
  for (unsigned int i = 0; i < _neighbourRanks.size(); i++) {
    std::cout << _neighbourRanks[i] << "  ";
  }
  std::cout << std::endl;
  std::cout << "Rank " << _rank << " has the following unique neighbours ";
  for (unsigned int i = 0; i < _numUniqueNeighbours; i++) {
    std::cout << _neighbourRanksUnique[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Rank " << _rank
            << " will be sending the following numbers of cells: ";
  for (unsigned int i = 0; i < _numUniqueNeighbours; i++) {
    std::cout << _numberOfCellsPerBuffer[i] << " ";
  }
  std::cout << std::endl;
#endif
}

void simplemd::services::ParallelTopologyService::initBuffers(
    const unsigned int &localNumberOfMolecules) {

  //idle processors should not enter this function, but still put a check:
  if (isIdle())
    return;

  unsigned int totalNumberOfLocalCells = 1;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    totalNumberOfLocalCells *= _localNumberOfCells[d];
  }

  // compute average number of molecules per linked cell
  double averageNumberOfMoleculesPerLinkedCell =
      (double) localNumberOfMolecules / (double) totalNumberOfLocalCells;

  // initialize ParallelAndLocalBufferService
  bool isOk;
  isOk =
      _bufferService.initialise(_numUniqueNeighbours, _numberOfCellsPerBuffer,
                                averageNumberOfMoleculesPerLinkedCell);
  if (!isOk) {
    std::cout
        << "Rank " << _rank
        << " could not initialise ParallelAndLocalBufferService. Terminating."
        << std::endl;
    exit(EXIT_FAILURE);
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "ParallelTopologyService: computed average number of molecules "
               "per linked cell: " << averageNumberOfMoleculesPerLinkedCell
            << std::endl;
#if (MD_PARALLEL == MD_YES)
  std::cout << "Rank " << _rank << " has the following buffer upper bounds: ";
  for (unsigned int i = 0; i < _numUniqueNeighbours; i++) {
    std::cout << _bufferService.getBufferCapacity(i) << " ";
  }
  std::cout << std::endl;
#endif
#endif

}

void simplemd::services::ParallelTopologyService::shutdown() {
  //All processors enter this function, but Idle processors do not have mpi
  //requests allocated
  // and their buffer services are not initialised, hence must not be shut down.
  if (!isIdle()) {
    //terminate Buffer services
    _bufferService.shutdown();
  }

  _globalIndexOfFirstCell.assign(0);
  _rank = 0;
}

int simplemd::services::ParallelTopologyService::getRank() const {
  return _rank;
}

const tarch::la::Vector<MD_DIM, unsigned int> &
simplemd::services::ParallelTopologyService::getProcessCoordinates() const {
  return _processCoordinates;
}

bool simplemd::services::ParallelTopologyService::isIdle() const {
  int processes = 1;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    processes = processes * _numberProcesses[d];
  }
  return (_rank > processes - 1);
}

std::vector<tarch::la::Vector<MD_DIM, unsigned int> >
simplemd::services::ParallelTopologyService::broadcastInnerCellViaBuffer(
    LinkedCell &cell, const unsigned int &cellIndex,
    const simplemd::services::LinkedCellService &linkedCellService) {

  std::vector<tarch::la::Vector<MD_DIM, unsigned int> > localIndex;

  // determine cell index in vector form
  tarch::la::Vector<MD_DIM, int> cellCoords(0);
  tarch::la::Vector<MD_DIM, unsigned int> size(
      linkedCellService.getLocalNumberOfCells() +
      2u * linkedCellService.getLocalIndexOfFirstCell());
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, tarch::la::Vector<MD_DIM, int> >
      neighbours(0);
  int help = (int) cellIndex;

#if (MD_DIM > 2)
  cellCoords[2] = help / (size[1] * size[0]);
  help = help - (cellCoords[2]) * (size[1] * size[0]);
#endif
#if (MD_DIM > 1)
  cellCoords[1] = help / (size[0]);
  help = help - cellCoords[1] * size[0];
#endif
  cellCoords[0] = help;
#if (MD_DEBUG == MD_YES)
  std::cout << "broadcastInnerCellViaBuffers: rank " << _rank
            << " cellCoords=" << cellCoords << std::endl;
#endif

#if (MD_DIM == 1)
  neighbours[0] = tarch::la::Vector<MD_DIM, int>(0);
  neighbours[1] = tarch::la::Vector<MD_DIM, int>(size(0) - 1);
#elif(MD_DIM == 2)
  neighbours[0] = tarch::la::Vector<MD_DIM, int>(0, 0);
  neighbours[1] = tarch::la::Vector<MD_DIM, int>(cellCoords[0], 0);
  neighbours[2] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, 0);
  neighbours[3] = tarch::la::Vector<MD_DIM, int>(0, cellCoords[1]);
  neighbours[4] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, cellCoords[1]);
  neighbours[5] = tarch::la::Vector<MD_DIM, int>(0, size[1] - 1);
  neighbours[6] = tarch::la::Vector<MD_DIM, int>(cellCoords[0], size[1] - 1);
  neighbours[7] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, size[1] - 1);
#elif(MD_DIM == 3)
  neighbours[0] = tarch::la::Vector<MD_DIM, int>(0, 0, 0);
  neighbours[1] = tarch::la::Vector<MD_DIM, int>(cellCoords[0], 0, 0);
  neighbours[2] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, 0, 0);

  neighbours[3] = tarch::la::Vector<MD_DIM, int>(0, cellCoords[1], 0);
  neighbours[4] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], cellCoords[1], 0);
  neighbours[5] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, cellCoords[1], 0);

  neighbours[6] = tarch::la::Vector<MD_DIM, int>(0, size[1] - 1, 0);
  neighbours[7] = tarch::la::Vector<MD_DIM, int>(cellCoords[0], size[1] - 1, 0);
  neighbours[8] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, size[1] - 1, 0);

  neighbours[9] = tarch::la::Vector<MD_DIM, int>(0, 0, cellCoords[2]);
  neighbours[10] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], 0, cellCoords[2]);
  neighbours[11] =
      tarch::la::Vector<MD_DIM, int>(size[0] - 1, 0, cellCoords[2]);

  neighbours[12] =
      tarch::la::Vector<MD_DIM, int>(0, cellCoords[1], cellCoords[2]);

  neighbours[13] =
      tarch::la::Vector<MD_DIM, int>(size[0] - 1, cellCoords[1], cellCoords[2]);

  neighbours[14] =
      tarch::la::Vector<MD_DIM, int>(0, size[1] - 1, cellCoords[2]);
  neighbours[15] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], size[1] - 1, cellCoords[2]);
  neighbours[16] =
      tarch::la::Vector<MD_DIM, int>(size[0] - 1, size[1] - 1, cellCoords[2]);

  neighbours[17] = tarch::la::Vector<MD_DIM, int>(0, 0, size[2] - 1);
  neighbours[18] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], 0, size[2] - 1);
  neighbours[19] = tarch::la::Vector<MD_DIM, int>(size[0] - 1, 0, size[2] - 1);

  neighbours[20] =
      tarch::la::Vector<MD_DIM, int>(0, cellCoords[1], size[2] - 1);
  neighbours[21] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], cellCoords[1], size[2] - 1);
  neighbours[22] =
      tarch::la::Vector<MD_DIM, int>(size[0] - 1, cellCoords[1], size[2] - 1);

  neighbours[23] = tarch::la::Vector<MD_DIM, int>(0, size[1] - 1, size[2] - 1);
  neighbours[24] =
      tarch::la::Vector<MD_DIM, int>(cellCoords[0], size[1] - 1, size[2] - 1);
  neighbours[25] =
      tarch::la::Vector<MD_DIM, int>(size[0] - 1, size[1] - 1, size[2] - 1);
#endif

  help = 0;
#if (MD_DIM > 2)
  for (int z = -1; z < 2; z++) {
#endif
#if (MD_DIM > 1)
    for (int y = -1; y < 2; y++) {
#endif
      for (int x = -1; x < 2; x++) {
        if (!((x == 0)
#if (MD_DIM > 1)
              && (y == 0)
#endif
#if (MD_DIM > 2)
              && (z == 0)
#endif
              )) {
          tarch::la::Vector<MD_DIM, int> diff(cellCoords - neighbours[help]);
          bool isNearBoundary = true;
          bool isInBoundary = true;
          for (unsigned int d = 0; d < MD_DIM; d++) {
            isNearBoundary = isNearBoundary && (abs(diff[d]) < 2);
            isInBoundary = isInBoundary && (diff[d] == 0);
          }
          // if this cell is close to a boundary and this boundary requires
          // communication...
          if (isNearBoundary && ((_boundary[help] == PERIODIC_BOUNDARY) ||
                                 (_boundary[help] == PARALLEL_BOUNDARY))) {
#if (MD_DEBUG == MD_YES)
            std::cout << "Cell is near boundary " << help << " ";
#endif
            int neighbourRank =
                ((x + ((int) _processCoordinates[0]) +
                  ((int) _numberProcesses[0])) % _numberProcesses[0])
#if (MD_DIM > 1)
                + ((y + ((int) _processCoordinates[1]) +
                    ((int) _numberProcesses[1])) % _numberProcesses[1]) *
                      _numberProcesses[0]
#endif
#if (MD_DIM > 2)
                + ((z + ((int) _processCoordinates[2]) +
                    ((int) _numberProcesses[2])) % _numberProcesses[2]) *
                      _numberProcesses[0] * _numberProcesses[1]
#endif
                ;
// if this is not the identity rank, do sending
#if (MD_DEBUG == MD_YES)
            std::cout << "NeighbourRank: " << neighbourRank
                      << " this rank: " << _rank << std::endl;
#endif
            if (neighbourRank != _rank) {
#if (MD_PARALLEL == MD_YES)
              // seek out current index for buffering
              int bufferIndex =
                  getCurrentBufferIndexFromNeighbourRank(neighbourRank);
#if (MD_DEBUG == MD_YES)
              std::cout << "Pack " << cell.getList().size()
                        << " molecules in buffer " << bufferIndex << std::endl;
#endif
              // push molecules into buffer
              for (std::list<Molecule *>::iterator it = cell.begin();
                   it != cell.end(); it++) {
                // for periodic boundaries: adapt position vector
                tarch::la::Vector<MD_DIM, double> position(
                    (*it)->getConstPosition());
                adaptPositionForPeriodicBoundaries(position, _boundary[help], x
#if (MD_DIM > 1)
                                                   ,
                                                   y
#endif
#if (MD_DIM > 2)
                                                   ,
                                                   z
#endif
                                                   );

#if (MD_DEBUG == MD_YES)
                std::cout << "Send molecule at position " << position
                          << ", velocity " << (*it)->getVelocity() << std::endl;
#endif
                // set correct position and push molecule
                pushMoleculeToSendBuffer(bufferIndex, *it, position);
              }
#else
              std::cout << "ERROR "
                           "simplemd::services::ParallelTopologyService::broadc"
                           "astInnerCell: There is no MPI support";
              std::cout << " available in the serial version ;-)" << std::endl;
              exit(EXIT_FAILURE);
#endif
              // if this cell is a ghost cell and periodic boundary condition
              // needs to be handled locally:
            } else if (neighbourRank == _rank && isInBoundary == true) {
#if (MD_DEBUG == MD_YES)
              std::cout << "Pack " << cell.getList().size()
                        << " molecules in local buffer " << std::endl;
#endif
              // push molecules into buffer
              for (std::list<Molecule *>::iterator it = cell.begin();
                   it != cell.end(); it++) {
                // for periodic boundaries: adapt position vector
                tarch::la::Vector<MD_DIM, double> position(
                    (*it)->getConstPosition());
                adaptPositionForPeriodicBoundaries(position, _boundary[help], x
#if (MD_DIM > 1)
                                                   ,
                                                   y
#endif
#if (MD_DIM > 2)
                                                   ,
                                                   z
#endif
                                                   );

#if (MD_DEBUG == MD_YES)
                std::cout << "Put molecule at position " << position
                          << ", velocity " << (*it)->getVelocity()
                          << ", forceOld" << (*it)->getForceOld()
                          << " in local buffer" << std::endl;
#endif
                // set correct position and push molecule
                pushMoleculeToLocalBuffer(*it, position);
              }
              // if this boundary is located on the same process:
            } else {
              // this is only supported for periodic boundaries...
              if (_boundary[help] != PERIODIC_BOUNDARY) {
                std::cout << "ERROR: "
                             "simplemd::services::ParallelTopologyService::broa"
                             "dcastInnerCellViaBuffers():";
                std::cout << " Parallel boundary detected where serial "
                             "boundary should be!" << std::endl;
                exit(EXIT_FAILURE);
              }
              // add local ghost cell coordinate where the molecules would need
              // to be copied to
              tarch::la::Vector<MD_DIM, unsigned int> serialNeighbour(0);
              for (unsigned int d = 0; d < MD_DIM; d++) {
                serialNeighbour[d] = (unsigned int)
                    neighbours[MD_LINKED_CELL_NEIGHBOURS - 1 - help][d];
              }
#if (MD_DEBUG == MD_YES)
              std::cout << "Serial neighbour found: Ghost cell coordinate: "
                        << serialNeighbour << std::endl;
#endif
              localIndex.push_back(serialNeighbour);
            }
          }

          // increment neighbor counter
          help++;
        }
      }
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif

  return localIndex;
}

bool simplemd::services::ParallelTopologyService::reduceGhostCellViaBuffer(
    LinkedCell &cell, const unsigned int &cellIndex,
    const simplemd::services::LinkedCellService &linkedCellService) {
#if (MD_PARALLEL == MD_YES)
  // determine neighbour rank (from position of the respective ghost cell
  // cellIndex)
  int cellCoords = 0;
  unsigned int help = cellIndex;
  int neighbourRank = 0;
#if (MD_DIM > 2)
  cellCoords = help / ((linkedCellService.getLocalNumberOfCells()[1] +
                        2 * linkedCellService.getLocalIndexOfFirstCell()[1]) *
                       (linkedCellService.getLocalNumberOfCells()[0] +
                        2 * linkedCellService.getLocalIndexOfFirstCell()[0]));
  help = help - ((unsigned int) cellCoords) *
                    ((linkedCellService.getLocalNumberOfCells()[1] +
                      2 * linkedCellService.getLocalIndexOfFirstCell()[1]) *
                     (linkedCellService.getLocalNumberOfCells()[0] +
                      2 * linkedCellService.getLocalIndexOfFirstCell()[0]));
  if (cellCoords == 0) {
    cellCoords = -1;
  } else if (cellCoords ==
             (int)(linkedCellService.getLocalNumberOfCells()[2] +
                   2 * linkedCellService.getLocalIndexOfFirstCell()[2] - 1)) {
    cellCoords = 1;
  } else {
    cellCoords = 0;
  }
  // the neighbouring process has coords cellCoords+_processCoordinates(2) (in
  // z-direction). To account for periodic
  // boundaries, we apply the modulo operator on top. As periodic boundaries at
  // the left/bottom,front would imply negative
  // coordinates that cannot be handled by the modulo operator correctly, we
  // further shift everything by _numberProcesses(2)
  // into the positive range.
  cellCoords = (((int)(_processCoordinates[2] + _numberProcesses[2])) +
                cellCoords) % _numberProcesses[2];
  neighbourRank += cellCoords * _numberProcesses[0] * _numberProcesses[1];
#endif
#if (MD_DIM > 1)
  cellCoords = help / (linkedCellService.getLocalNumberOfCells()[0] +
                       2 * linkedCellService.getLocalIndexOfFirstCell()[0]);
  help =
      help - cellCoords * (linkedCellService.getLocalNumberOfCells()[0] +
                           2 * linkedCellService.getLocalIndexOfFirstCell()[0]);
  if (cellCoords == 0) {
    cellCoords = -1;
  } else if (cellCoords ==
             (int)(linkedCellService.getLocalNumberOfCells()[1] +
                   2 * linkedCellService.getLocalIndexOfFirstCell()[1] - 1)) {
    cellCoords = 1;
  } else {
    cellCoords = 0;
  }
  cellCoords = (((int)(_processCoordinates[1] + _numberProcesses[1])) +
                cellCoords) % _numberProcesses[1];
  neighbourRank += cellCoords * _numberProcesses[0];
#endif
  cellCoords = help;
  if (cellCoords == 0) {
    cellCoords = -1;
  } else if (cellCoords ==
             (int)(linkedCellService.getLocalNumberOfCells()[0] +
                   2 * linkedCellService.getLocalIndexOfFirstCell()[0] - 1)) {
    cellCoords = 1;
  } else {
    cellCoords = 0;
  }
  cellCoords = (((int)(_processCoordinates[0] + _numberProcesses[0])) +
                cellCoords) % _numberProcesses[0];
  neighbourRank += cellCoords;

  // if this is the same process, return false
  if (neighbourRank == _rank) {
    return false;
  }

  // if this is a neighbour, push in send buffer
  if (isParallelNeighbour(neighbourRank)) {

    // determine appropriate buffer index
    int bufferIndex = getCurrentBufferIndexFromNeighbourRank(neighbourRank);
    for (std::list<Molecule *>::iterator it = cell.begin(); it != cell.end();
         it++) {
#if (MD_DEBUG == MD_YES)
      std::cout << "Reduce molecule at position " << (*it)->getConstPosition()
                << ", velocity " << (*it)->getConstVelocity() << std::endl;
#endif
      pushMoleculeToSendBuffer(bufferIndex, *it, (*it)->getConstPosition());
    }
    return true;
  } else {
    return false;
  }
#else
  return false;
#endif
}

void simplemd::services::ParallelTopologyService::unpackLocalBuffer(
    simplemd::services::MoleculeService &moleculeService,
    simplemd::services::LinkedCellService &linkedCellService) {
  unpackBuffer(_bufferService.getLocalBuffer(), moleculeService,
               linkedCellService);
}

void simplemd::services::ParallelTopologyService::communicationSteps_1_2() {
#if (MD_PARALLEL == MD_YES)
  //issue all Irecv calls
  for (unsigned int bufferIndex = 0; bufferIndex < _numUniqueNeighbours;
       bufferIndex++) {
    bufferIrecv(_bufferService.getReceiveBuffer(bufferIndex),
                _neighbourRanksUnique[bufferIndex],
                _receiveRequests[bufferIndex]);
  }

  //issue all Isend calls
  for (unsigned int bufferIndex = 0; bufferIndex < _numUniqueNeighbours;
       bufferIndex++) {
    bufferIsend(_bufferService.getSendBuffer(bufferIndex),
                _neighbourRanksUnique[bufferIndex], _sendRequests[bufferIndex]);
  }
#endif
}

void simplemd::services::ParallelTopologyService::communicationSteps_3_4(
    simplemd::services::MoleculeService &moleculeService,
    simplemd::services::LinkedCellService &linkedCellService) {
#if (MD_PARALLEL == MD_YES)
  // Steps 3 and 4 together:

  MPI_Status stat;
  int result;
  int flag = 0;
  int count = 0;
  unsigned int i_buf;

  // flags to indicate which request has been retired
  bool allCompleted;
  bool recvCompleted[MD_LINKED_CELL_NEIGHBOURS];
  bool sendCompleted[MD_LINKED_CELL_NEIGHBOURS];

  // set flags accordingly
  for (i_buf = 0; i_buf < _numUniqueNeighbours; i_buf++) {
    recvCompleted[i_buf] = false;
    sendCompleted[i_buf] = false;
  }
  allCompleted = false;

  while (!allCompleted) {
    allCompleted = true;

    for (i_buf = 0; i_buf < _numUniqueNeighbours; i_buf++) {
      //test send request
      if (!sendCompleted[i_buf]) {
        result = MPI_Test(&(_sendRequests[i_buf]), &flag, &stat);
        if (result != MPI_SUCCESS) {
          std::cout << "rank " << _rank << " was not able to test send request "
                    << i_buf << std::endl;
          exit(EXIT_FAILURE);
        }
        // if flag == 1, this means message request has just been completed
        if (flag == 1) {
          //clear buffer
          (_bufferService.getSendBuffer(i_buf))->clearBuffer();
          sendCompleted[i_buf] = true;
#if (MD_DEBUG == MD_YES)
          std::cout << "Rank " << _rank << " sent "
                    << _bufferService.getSendBufferLength(i_buf) / MD_DIM / 3
                    << " molecules to " << _neighbourRanksUnique[i_buf]
                    << std::endl;
#endif
        }
      }
      allCompleted = allCompleted && sendCompleted[i_buf];

      //test receive request
      if (!recvCompleted[i_buf]) {
        result = MPI_Test(&(_receiveRequests[i_buf]), &flag, &stat);
        if (result != MPI_SUCCESS) {
          std::cout << "rank " << _rank
                    << " was not able to test receive request " << i_buf
                    << std::endl;
          exit(EXIT_FAILURE);
        }
        // if flag == 1, this means message request has just been completed
        if (flag == 1) {
          //set length of receive buffer correctly
          MPI_Get_count(&stat, MPI_DOUBLE, &count);
          _bufferService.setReceiveBufferLength(i_buf, (unsigned int) count);
          // DO NOT unpack buffer yet - this would indeed be an additional
          // overlap of waiting time,
          // but it makes the unpacking order random, which was observed to
          // affect (change) roundoffs between subsequent runs for a fixed
          // randomization seedpoint.
          // Unpacking all buffers after reception has been completed gives a
          // well defined insertion
          // order and makes parallel runs comparable to each other

          // An alternative implementation is present in older revisions of the
          // repository
          recvCompleted[i_buf] = true;
#if (MD_DEBUG == MD_YES)
          std::cout << "Rank " << _rank << " received " << count / MD_DIM / 3
                    << " molecules from " << _neighbourRanksUnique[i_buf]
                    << std::endl;
#endif
        }
      }
      allCompleted = allCompleted && recvCompleted[i_buf];

    }

  } // end of while loop

  // now unpack all buffers in a fixed order:
  for (i_buf = 0; i_buf < _numUniqueNeighbours; i_buf++) {
    unpackBuffer(_bufferService.getReceiveBuffer(i_buf), moleculeService,
                 linkedCellService);
  }

#endif
}

bool simplemd::services::ParallelTopologyService::globalToLocalRegionOfInterest(
    const tarch::la::Vector<MD_DIM, unsigned int> &globalStartCell,
    const tarch::la::Vector<MD_DIM, unsigned int> &globalRange,
    tarch::la::Vector<MD_DIM, unsigned int> &localStartCell,
    tarch::la::Vector<MD_DIM, unsigned int> &localRange) const {

  bool isIntersecting = true;
  int help; //need an int because some entries will be negative

  // for convenience, we will work with a start point and an endpoint
  tarch::la::Vector<MD_DIM, unsigned int> localEndCell;
  tarch::la::Vector<MD_DIM, unsigned int> globalEndCell;

  for (unsigned int d = 0; d < MD_DIM; d++) {
    if (globalRange[d] == 0) {
      isIntersecting = false;
      break;
    }
    globalEndCell[d] = globalStartCell[d] + globalRange[d] - 1;

    //compute start point

    help = (int) globalStartCell[d] -
           (int)(_processCoordinates[d] * _localNumberOfCells[d]);
    // first, set negative entries to 0
    help = std::max(help, 0);
    // now if entry is 0, but process is not on boundary,
    // then this cell in fact belongs to the previous process along the same
    // dimension,
    // so set to 1. (In order for regions to be non-overlapping)
    if (help == 0 && _processCoordinates[d] != 0)
      help = 1;
    // this is now the correct value in the local frame of reference.
    // now if the starting point is not inside the domain, there is no
    // intersection
    if ((unsigned int) help > _localNumberOfCells[d]) {
      isIntersecting = false;
      break;
    }
    localStartCell[d] = (unsigned int) help;

    //compute end point
    help = (int) globalEndCell[d] -
           (int)(_processCoordinates[d] * _localNumberOfCells[d]);
    // first, set values larger than the number of cells along the dimension +1
    // to the number of cells +1 (+1 is needed because of ghost layer)
    help = std::min(help, (int)(_localNumberOfCells[d] + 1));
    // now if entry is equal to _localNumberOfCells(d) + 1, but process is not
    // on boundary,
    // then this cell in fact belongs to the next process along the same
    // dimension
    // so set to _localNumberOfCells(d). (in order for regions to be
    // non-overlapping)
    if ((help == (int)(_localNumberOfCells[d] + 1)) &&
        _processCoordinates[d] != (_numberProcesses[d] - 1))
      help = _localNumberOfCells[d];
    // this is now the correct value in the local frame of reference.
    // now if the end point is not inside the domain, there is no intersection
    if (help < 1) {
      isIntersecting = false;
      break;
    }
    localEndCell[d] = (unsigned int) help;

    //compute local range
    localRange[d] = localEndCell[d] - localStartCell[d] + 1;
  }

  // at this point, if there is an intersection, everything is already
  // calculated
  if (isIntersecting)
    return true;
  else {
    // if not, set output vectors to 0
    localStartCell = tarch::la::Vector<MD_DIM, unsigned int>(0);
    localRange = tarch::la::Vector<MD_DIM, unsigned int>(0);
    return false;
  }
}

tarch::la::Vector<MD_DIM, unsigned int>
simplemd::services::ParallelTopologyService::localToGlobalCellIndexVector(
    const tarch::la::Vector<MD_DIM, unsigned int> &localCellIndexVector) const {
  tarch::la::Vector<MD_DIM, unsigned int> globalCellIndexVector;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    globalCellIndexVector[d] = localCellIndexVector[d] +
                               _processCoordinates[d] * _localNumberOfCells[d];
  }
  return globalCellIndexVector;
}

void simplemd::services::ParallelTopologyService::createNeighbourRanks(
    const tarch::la::Vector<MD_DIM, unsigned int> &processCoordinates,
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &
        localBoundary,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
    std::vector<int> &neighbourRanks, std::vector<int> &neighbourRanksUnique,
    unsigned int &numUniqueNeighbours,
    unsigned int numberOfCellsPerBuffer[]) const {
  unsigned int neighbourCounter = 0;
  tarch::la::Vector<MD_DIM, int> tmpCoords(0);
  // empty ranks vector and reset number of transferred linked cells
  neighbourRanks.clear();
  neighbourRanksUnique.clear();
  numUniqueNeighbours = 0;

  //set initial upper bounds to zero
  for (unsigned int i = 0; i < MD_LINKED_CELL_NEIGHBOURS; i++) {
    numberOfCellsPerBuffer[i] = 0;
  }

#if (MD_DIM > 2)
  for (int z = -1; z < 2; z++) {
    tmpCoords[2] = ((int) processCoordinates[2]) + z;
#endif
#if (MD_DIM > 1)
    for (int y = -1; y < 2; y++) {
      tmpCoords[1] = ((int) processCoordinates[1]) + y;
#endif
      for (int x = -1; x < 2; x++) {
        tmpCoords[0] = ((int) processCoordinates[0]) + x;
        if (!((x == 0)
#if (MD_DIM > 1)
              && (y == 0)
#endif
#if (MD_DIM > 2)
              && (z == 0)
#endif
              )) {
          int tmpRank = 0;
          bool isBoundary = false;
          // find out if this process would be outside the domain and correct
          // its position in case of
          // periodic boundaries
          for (unsigned int d = 0; d < MD_DIM; d++) {
            if (tmpCoords[d] < 0) {
              tmpCoords[d] += numberProcesses[d];
              isBoundary = true;
            } else if (tmpCoords[d] > (int)(numberProcesses[d] - 1)) {
              tmpCoords[d] -= numberProcesses[d];
              isBoundary = true;
            }
          }

          // if this neighbour corresponds to a periodic-boundary-process or an
          // inner process,
          // add the rank and the respective number of transferred cells
          if ((isBoundary && (localBoundary[neighbourCounter] ==
                              simplemd::PERIODIC_BOUNDARY)) || (!isBoundary)) {
            // compute and add rank
            tmpRank = tmpCoords[0]
#if (MD_DIM > 1)
                      + tmpCoords[1] * numberProcesses[0]
#endif
#if (MD_DIM > 2)
                      + tmpCoords[2] * numberProcesses[1] * numberProcesses[0]
#endif
                ;

            // if this is not the same process, add respective information
            if (tmpRank != _rank) {
              neighbourRanks.push_back(tmpRank);
              addNeighbourToNeighbourRanksUnique(neighbourRanksUnique,
                                                 numUniqueNeighbours, tmpRank);

              int bufIndex = getCurrentBufferIndexFromNeighbourRank(tmpRank);
              // compute and add cells to be transferred from this neighbour
              int numCells = getNumberOfTransferredCells(x
#if (MD_DIM > 1)
                                                         ,
                                                         y
#endif
#if (MD_DIM > 2)
                                                         ,
                                                         z
#endif
                                                         );
              /*the following value should not be needed anymore, to be removed
               * soon*/
              numberOfCellsPerBuffer[bufIndex] += (unsigned int) numCells;

            }
          }

          // increment neighbour counter
          neighbourCounter++;

          //undo effects of loop that sets isBoundary -
          // otherwise we proceed subsequent checks with possibly modified
          // tmpCoords
          // which leads to incorrect setting of isBoundary flag
          tmpCoords[0] = ((int) processCoordinates[0]) + x;
#if (MD_DIM > 1)
          tmpCoords[1] = ((int) processCoordinates[1]) + y;
#endif
#if (MD_DIM > 2)
          tmpCoords[2] = ((int) processCoordinates[2]) + z;
#endif
        }
      }
#if (MD_DIM > 1)
    }
#endif
#if (MD_DIM > 2)
  }
#endif
}

// add ranks one by one, skipping repeating ones
void
simplemd::services::ParallelTopologyService::addNeighbourToNeighbourRanksUnique(
    std::vector<int> &neighbourRanksUnique, unsigned int &numUniqueNeighbours,
    const int &addedNeighbour) const {
  //note that due to functionality of createNeighbourRanks, _rank cannot appear
  //in this vector
  bool repetitionFound = false;
  for (unsigned int j = 0; j < numUniqueNeighbours; j++) {
    if (addedNeighbour == neighbourRanksUnique[j]) {
      repetitionFound = true;
      break;
    }
  }
  if (!repetitionFound) {
    neighbourRanksUnique.push_back(addedNeighbour);
    numUniqueNeighbours++;
  }
}

tarch::la::Vector<MD_DIM, double>
simplemd::services::ParallelTopologyService::computeMeshwidth(
    const tarch::la::Vector<MD_DIM, double> &prescribedWidth,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
    const tarch::la::Vector<MD_DIM, double> &domainSize) const {
  tarch::la::Vector<MD_DIM, double> meshWidth;
  for (unsigned int d = 0; d < MD_DIM; d++) {
    // compute approx. number of cells in d-direction
    unsigned int cells = static_cast<unsigned int>(
        floor((domainSize[d] / numberProcesses[d]) / prescribedWidth[d] + 0.5));
    // determine correct mesh width in d-direction when using "cells" cells
    meshWidth[d] = (domainSize[d] / numberProcesses[d]) / cells;
  }
#if (MD_DEBUG == MD_YES)
  std::cout << "Meshwidth: " << meshWidth << std::endl;
#endif
  return meshWidth;
}

tarch::la::Vector<MD_DIM, unsigned int>
simplemd::services::ParallelTopologyService::computeNumberOfCells(
    const tarch::la::Vector<MD_DIM, double> &meshWidth,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
    const tarch::la::Vector<MD_DIM, double> &domainSize) const {
  tarch::la::Vector<MD_DIM, unsigned int> numberOfCells(0);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    // compute approx. number of cells in d-direction
    unsigned int cells = static_cast<unsigned int>(
        floor((domainSize[d] / numberProcesses[d]) / meshWidth[d] + 0.5));
    // determine correct mesh width in d-direction when using "cells" cells
    double tmpWidth = (domainSize[d] / numberProcesses[d]) / cells;
    // compute real number of cells
    numberOfCells[d] = static_cast<unsigned int>(
        floor((domainSize[d] / numberProcesses[d]) / tmpWidth + 0.5));
  }
#if (MD_DEBUG == MD_YES)
  std::cout << "Number cells: " << numberOfCells << std::endl;
#endif
  return numberOfCells;
}

tarch::la::Vector<MD_DIM, unsigned int>
simplemd::services::ParallelTopologyService::computeGlobalNumberOfCells(
    const tarch::la::Vector<MD_DIM, double> &meshWidth,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses,
    const tarch::la::Vector<MD_DIM, double> &domainSize) const {
  tarch::la::Vector<MD_DIM, unsigned int> globalCells =
      computeNumberOfCells(meshWidth, numberProcesses, domainSize);
  for (unsigned int d = 0; d < MD_DIM; d++) {
    globalCells[d] *= numberProcesses[d];
  }
  return globalCells;
}

bool simplemd::services::ParallelTopologyService::isParallelNeighbour(
    const int &neighbourRank) const {
  for (unsigned int i = 0; i < _numUniqueNeighbours; i++) {
    if (neighbourRank == _neighbourRanksUnique[i]) {
      return true;
    }
  }
  return false;
}

tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>
simplemd::services::ParallelTopologyService::computeLocalBoundaryInformation(
    const tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType> &
        boundary,
    const tarch::la::Vector<MD_DIM, unsigned int> &processCoordinates,
    const tarch::la::Vector<MD_DIM, unsigned int> &numberProcesses) {
  tarch::la::Vector<MD_LINKED_CELL_NEIGHBOURS, simplemd::BoundaryType>
      localBoundary(PARALLEL_BOUNDARY);
  unsigned int counter = 0;
#if (MD_DIM == 1)
  for (int x = ((int) processCoordinates[0]) - 1;
       x < ((int) processCoordinates[0]) + 2; x++) {
    if (!(x == (int) processCoordinates[0])) {
      tarch::la::Vector<MD_DIM, int> nextProcess(x);
      // left boundary
      if (nextProcess[0] == -1) {
        localBoundary[counter] = boundary[0];
        // right boundary
      } else if (nextProcess[0] == (int) numberProcesses[0]) {
        localBoundary[counter] = boundary[1];
        // parallel boundaries in the middle of domain
      } else {
        localBoundary[counter] = simplemd::PARALLEL_BOUNDARY;
      }

      counter++;
    }
  }
#elif(MD_DIM == 2)
  for (int y = ((int) processCoordinates[1]) - 1;
       y < ((int) processCoordinates[1]) + 2; y++) {
    for (int x = ((int) processCoordinates[0]) - 1;
         x < ((int) processCoordinates[0]) + 2; x++) {
      if (!((y == (int) processCoordinates[1]) &&
            (x == (int) processCoordinates[0]))) {
        tarch::la::Vector<MD_DIM, int> nextProcess(x, y);

        // lower left corner
        if ((nextProcess[0] == -1) && (nextProcess[1] == -1)) {
          localBoundary[counter] = boundary[0];
          // lower right corner
        } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                   (nextProcess[1] == -1)) {
          localBoundary[counter] = boundary[2];
          // upper left corner
        } else if ((nextProcess[0] == -1) &&
                   (nextProcess[1] == (int) numberProcesses[1])) {
          localBoundary[counter] = boundary[5];
        } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                   (nextProcess[1] == (int) numberProcesses[1])) {
          localBoundary[counter] = boundary[7];
          // left edge
        } else if (nextProcess[0] == -1) {
          localBoundary[counter] = boundary[3];
          // right edge
        } else if (nextProcess[0] == (int) numberProcesses[0]) {
          localBoundary[counter] = boundary[4];
          // lower edge
        } else if (nextProcess[1] == -1) {
          localBoundary[counter] = boundary[1];
          // upper edge
        } else if (nextProcess[1] == (int) numberProcesses[1]) {
          localBoundary[counter] = boundary[6];
          // mid boundary
        } else {
          localBoundary[counter] = simplemd::PARALLEL_BOUNDARY;
        }

        counter++;
      }
    }
  }
#elif(MD_DIM == 3)
  for (int z = ((int) processCoordinates[2]) - 1;
       z < ((int) processCoordinates[2]) + 2; z++) {
    for (int y = ((int) processCoordinates[1]) - 1;
         y < ((int) processCoordinates[1]) + 2; y++) {
      for (int x = ((int) processCoordinates[0]) - 1;
           x < ((int) processCoordinates[0]) + 2; x++) {
        if (!((z == (int) processCoordinates[2]) &&
              (y == (int) processCoordinates[1]) &&
              (x == (int) processCoordinates[0]))) {
          tarch::la::Vector<MD_DIM, int> nextProcess(x, y, z);

          // -------------- corners -------------------------

          // lower left front corner
          if ((nextProcess[0] == -1) && (nextProcess[1] == -1) &&
              (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[0];
            // lower right front corner
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == -1) && (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[2];
            // lower left back corner
          } else if ((nextProcess[0] == -1) &&
                     (nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[6];
            // lower right back corner
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[8];
            // upper left front corner
          } else if ((nextProcess[0] == -1) && (nextProcess[1] == -1) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[17];
            // upper right front corner
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == -1) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[19];
            // upper left back corner
          } else if ((nextProcess[0] == -1) &&
                     (nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[23];
            // upper right back corner
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[25];
          }

              // --------------------- edges --------------
              // ----------- z-aligned ------------
              // left front edge
              else if ((nextProcess[0] == -1) && (nextProcess[1] == -1)) {
            localBoundary[counter] = boundary[9];
            // right front edge
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == -1)) {
            localBoundary[counter] = boundary[11];
            // left back edge
          } else if ((nextProcess[0] == -1) &&
                     (nextProcess[1] == (int) numberProcesses[1])) {
            localBoundary[counter] = boundary[14];
            // right back edge
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[1] == (int) numberProcesses[1])) {
            localBoundary[counter] = boundary[16];
          }
              // ----------- y-aligned ------------
              // lower left edge
              else if ((nextProcess[0] == -1) && (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[3];
            // lower right edge
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[5];
            // upper left edge
          } else if ((nextProcess[0] == -1) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[20];
            // upper right edge
          } else if ((nextProcess[0] == (int) numberProcesses[0]) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[22];
          }
              // ----------- x-aligned ------------
              // lower front edge
              else if ((nextProcess[1] == -1) && (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[1];
            // lower back edge
          } else if ((nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == -1)) {
            localBoundary[counter] = boundary[7];
            // upper front edge
          } else if ((nextProcess[1] == -1) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[18];
            // upper back edge
          } else if ((nextProcess[1] == (int) numberProcesses[1]) &&
                     (nextProcess[2] == (int) numberProcesses[2])) {
            localBoundary[counter] = boundary[24];
          }

              // ------------------- faces ---------------
              // left face
              else if (nextProcess[0] == -1) {
            localBoundary[counter] = boundary[12];
            // right face
          } else if (nextProcess[0] == (int) numberProcesses[0]) {
            localBoundary[counter] = boundary[13];
            // front face
          } else if (nextProcess[1] == -1) {
            localBoundary[counter] = boundary[10];
            // back face
          } else if (nextProcess[1] == (int) numberProcesses[1]) {
            localBoundary[counter] = boundary[15];
            // bottom face
          } else if (nextProcess[2] == -1) {
            localBoundary[counter] = boundary[4];
            // top face
          } else if (nextProcess[2] == (int) numberProcesses[2]) {
            localBoundary[counter] = boundary[21];
            // otherwise: parallel boundary
          } else {
            localBoundary[counter] = simplemd::PARALLEL_BOUNDARY;
          }

          counter++;
        }
      } // x
    }   // y
  }     // z
#endif
  return localBoundary;
}

unsigned int
simplemd::services::ParallelTopologyService::getNumberOfTransferredCells(
    const int &x
#if (MD_DIM > 1)
    ,
    const int &y
#endif
#if (MD_DIM > 2)
    ,
    const int &z
#endif
    ) const {
  unsigned int cells = 0;
#if (MD_DIM == 1)
#if (MD_ERROR == MD_YES)
  if ((x != 1) && (x != -1)) {
    std::cout << "ERROR "
                 "simplemd::services::ParallelTopologyService::getNumberOfTrans"
                 "ferredCells: x != 1 and x != -1" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
  return 1;
#elif(MD_DIM == 2)
#if (MD_ERROR == MD_YES)
  int sum = abs(x) + abs(y);
  if ((sum != 1) && (sum != 2)) {
    std::cout << "ERROR "
                 "simplemd::services::ParallelTopologyService::getNumberOfTrans"
                 "ferredCells: Unknown neighbour relation!" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif
  // left/ right edge
  if (((x == -1) && (y == 0)) || ((x == 1) && (y == 0))) {
    return _localNumberOfCells[1];
    // front/ back edge
  } else if (((y == -1) && (x == 0)) || ((y == 1) && (x == 0))) {
    return _localNumberOfCells[0];
    // otherwise, this is a corner
  } else {
    return 1;
  }
#elif(MD_DIM == 3)
  // corners
  if (abs(x) + abs(y) + abs(z) == 3) {
    return 1;
    // x-axis aligned edges
  } else if ((x == 0) && (((y == -1) && (z == -1)) || ((y == -1) && (z == 1)) ||
                          ((y == 1) && (z == -1)) || ((y == 1) && (z == 1)))) {
    return _localNumberOfCells[0];
    // y-axis aligned edges
  } else if ((y == 0) && (((x == -1) && (z == -1)) || ((x == -1) && (z == 1)) ||
                          ((x == 1) && (z == -1)) || ((x == 1) && (z == 1)))) {
    return _localNumberOfCells[1];
    // z-axis aligned edges
  } else if ((z == 0) && (((y == -1) && (x == -1)) || ((y == -1) && (x == 1)) ||
                          ((y == 1) && (x == -1)) || ((y == 1) && (x == 1)))) {
    return _localNumberOfCells[2];
    // left/ right plane
  } else if ((y == 0) && (z == 0) && ((x == -1) || (x == 1))) {
    return _localNumberOfCells[1] * _localNumberOfCells[2];
    // front/back plane
  } else if ((x == 0) && (z == 0) && ((y == -1) || (y == 1))) {
    return _localNumberOfCells[0] * _localNumberOfCells[2];
    // top/ bottom plane
  } else if ((x == 0) && (y == 0) && ((z == -1) || (z == 1))) {
    return _localNumberOfCells[0] * _localNumberOfCells[1];
  } else {
#if (MD_ERROR == MD_YES)
    std::cout << "ERROR "
                 "simplemd::services::ParallelTopologyService::getNumberOfTrans"
                 "ferredCells: Unknown neighbour relation!" << std::endl;
    exit(EXIT_FAILURE);
#endif
    return 0;
  }
#endif
  return cells;
}

void
simplemd::services::ParallelTopologyService::adaptPositionForPeriodicBoundaries(
    tarch::la::Vector<MD_DIM, double> &position,
    const simplemd::BoundaryType &boundaryType, const int &x
#if (MD_DIM > 1)
    ,
    const int &y
#endif
#if (MD_DIM > 2)
    ,
    const int &z
#endif
    ) const {

  if (boundaryType == PERIODIC_BOUNDARY) {
    if (((int) _processCoordinates[0]) + x < 0) {
      position[0] += _domainSize[0];
    }
    if (((int) _processCoordinates[0]) + x > ((int) _numberProcesses[0] - 1)) {
      position[0] -= _domainSize[0];
    }

#if (MD_DIM > 1)
    if (((int) _processCoordinates[1]) + y < 0) {
      position[1] += _domainSize[1];
    }
    if (((int) _processCoordinates[1]) + y > ((int) _numberProcesses[1] - 1)) {
      position[1] -= _domainSize[1];
    }
#endif

#if (MD_DIM > 2)
    if (((int) _processCoordinates[2]) + z < 0) {
      position[2] += _domainSize[2];
    }
    if (((int) _processCoordinates[2]) + z > ((int) _numberProcesses[2] - 1)) {
      position[2] -= _domainSize[2];
    }
#endif
  }
}

// get current index of buffer arrays from neighbourRank
unsigned int simplemd::services::ParallelTopologyService::
    getCurrentBufferIndexFromNeighbourRank(const int &neighbourRank) const {
  unsigned int bufferIndex;
  bool found = false;
  for (bufferIndex = 0; bufferIndex < _numUniqueNeighbours; bufferIndex++) {
    if (neighbourRank == _neighbourRanksUnique[bufferIndex]) {
      found = true;
      break;
    }
  }

  if (!found) {
    std::cout << "Problems: "
                 "ParallelTopologyService.cpp::getCurrentBufferIndexFromNeighbo"
                 "urRank() no corresponding index found" << std::endl;
    exit(EXIT_FAILURE);
    return 0;
  } else {
    return bufferIndex;
  }
}

void simplemd::services::ParallelTopologyService::unpackBuffer(
    ParallelAndLocalBufferService::SimpleBuffer *buf,
    simplemd::services::MoleculeService &moleculeService,
    simplemd::services::LinkedCellService &linkedCellService) {
  tarch::la::Vector<MD_DIM, double> position(0.0);
  tarch::la::Vector<MD_DIM, double> velocity(0.0);
  tarch::la::Vector<MD_DIM, double> forceOld(0.0);
  tarch::la::Vector<MD_DIM, unsigned int> cellIndex(0);
  bool isFixed;

  unsigned int iMol;

  const double *const values = buf->getValues();
  const unsigned int len = buf->getLength();

  for (iMol = 0; iMol < len; iMol += (3 * MD_DIM + 1)) {
// unroll loops manually:
#if (MD_DIM == 1)
    position[0] = values[iMol];
    velocity[0] = values[iMol + 1];
    forceOld[0] = values[iMol + 2];
    isFixed = (bool) values[iMol + 3];
#endif
#if (MD_DIM == 2)
    position[0] = values[iMol];
    position[1] = values[iMol + 1];
    velocity[0] = values[iMol + 2];
    velocity[1] = values[iMol + 3];
    forceOld[0] = values[iMol + 4];
    forceOld[1] = values[iMol + 5];
    isFixed = (bool) values[iMol + 6];
#endif
#if (MD_DIM == 3)
    position[0] = values[iMol];
    position[1] = values[iMol + 1];
    position[2] = values[iMol + 2];
    velocity[0] = values[iMol + 3];
    velocity[1] = values[iMol + 4];
    velocity[2] = values[iMol + 5];
    forceOld[0] = values[iMol + 6];
    forceOld[1] = values[iMol + 7];
    forceOld[2] = values[iMol + 8];
    isFixed = (bool) values[iMol + 9];
#endif

    // determine cell index vector
    for (unsigned int d = 0; d < MD_DIM; d++) {
      int cell = (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d]));
      cell -= (int)(_processCoordinates[d] * _localNumberOfCells[d]);
      cell += (int) linkedCellService.getLocalIndexOfFirstCell()[d];
#if (MD_ERROR == MD_YES)
      if (cell < 0) {
        std::cout
            << "ERROR "
               "simplemd::services::ParallelTopologyService::unpackBuffer:";
        std::cout << " Molecule out of range: Position=" << position
                  << ", cell=" << cell;
        std::cout << " _processCoords=" << _processCoordinates
                  << ", loc.num=" << _localNumberOfCells
                  << ", domOff=" << _domainOffset;
        std::cout << ", mesh=" << _meshWidth << std::endl;
        exit(EXIT_FAILURE);
      }
#endif
      cellIndex[d] = (unsigned int) cell;
    }

    // add molecule to MoleculeService and LinkedCellService
    Molecule myMolecule(position, velocity);
    myMolecule.setForceOld(forceOld);
    if (isFixed)
      myMolecule.fix();

    Molecule *mPtr = moleculeService.addMolecule(myMolecule);
#if (MD_DEBUG == MD_YES)
    std::cout << "Rank " << _rank
              << ": unpacked molecule from buffer into cell " << cellIndex
              << std::endl;
#endif
    linkedCellService.addMoleculeToLinkedCell(*mPtr, cellIndex);
  }

  // clear the buffer
  buf->clearBuffer();

}

void simplemd::services::ParallelTopologyService::pushMoleculeToLocalBuffer(
    const Molecule *mol, const tarch::la::Vector<MD_DIM, double> &pos) {
#if (MD_ERROR == MD_YES)
  bool isOk;
  isOk =
#endif
      _bufferService.pushMoleculeToLocalBuffer(pos, mol);
#if (MD_ERROR == MD_YES)
  //if errors are off, application will terminate earlier
  if (!isOk) {
    std::cout << "Rank " << _rank
              << " was unable to push molecule to local buffer. " << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

#if (MD_DEBUG == MD_YES)
  std::cout
      << "rank " << _rank
      << " pushed molecule to local buffer successfully, buffer length is "
      << _bufferService.getLocalBuffer()->getLength() << std::endl;
  std::cout << "pushed molecule with position " << pos << ", velocity "
            << mol->getConstVelocity() << ", forceOld"
            << mol->getConstForceOld() << std::endl;
#endif
}

#if (MD_PARALLEL == MD_YES)
void simplemd::services::ParallelTopologyService::pushMoleculeToSendBuffer(
    const unsigned int &bufferIndex, const Molecule *mol,
    const tarch::la::Vector<MD_DIM, double> &pos) {
#if (MD_ERROR == MD_YES)
  bool isOk;
  isOk =
#endif
      _bufferService.pushMoleculeToSendBuffer(pos, mol, bufferIndex);
#if (MD_ERROR == MD_YES)
  //if errors are off, application will terminate earlier
  if (!isOk) {
    std::cout << "Rank " << _rank
              << " was unable to push molecule to send buffer " << bufferIndex
              << ". " << std::endl;
    std::cout << "Upper bound of buffer ("
              << _bufferService.getBufferCapacity(bufferIndex)
              << ") exceeded. Fatal error. Please contact developers as this "
                 "is a principal issue." << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

#if (MD_DEBUG == MD_YES)
  std::cout << "Rank " << _rank << " pushed molecule to send buffer "
            << bufferIndex << " successfully, buffer length is "
            << _bufferService.getSendBufferLength(bufferIndex) << std::endl;
  std::cout << "pushed molecule at position " << pos << ", velocity "
            << mol->getConstVelocity() << ", forceOld"
            << mol->getConstForceOld() << std::endl;
#endif
}
#endif

#if (MD_PARALLEL == MD_YES)
void simplemd::services::ParallelTopologyService::bufferIsend(
    ParallelAndLocalBufferService::SimpleBuffer *buffer,
    const int &neighbourRank, MPI_Request &request) const {
  int result;

  // send only as much values as we need to
  result = MPI_Isend(buffer->getValues(), buffer->getLength(), MPI_DOUBLE,
                     neighbourRank, _bufferTag, _communicator, &request);

  if (result != MPI_SUCCESS) {
    std::cout << "rank " << _rank << " was not able to Isend buffer to "
              << neighbourRank << ": " << std::endl;
    //std::cout << tarch::parallel::MPIReturnValueToString(result);
    exit(EXIT_FAILURE);
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "Rank " << _rank << " is sending "
            << buffer->getLength() / 3 / MD_DIM << " molecules to "
            << neighbourRank << std::endl;
#endif
}
#endif

#if (MD_PARALLEL == MD_YES)
void simplemd::services::ParallelTopologyService::bufferIrecv(
    ParallelAndLocalBufferService::SimpleBuffer *buffer,
    const int &neighbourRank, MPI_Request &request) const {
  int result;

  // receive as many molecules as the upper bound is
  result = MPI_Irecv(buffer->getValues(), buffer->getCapacity(), MPI_DOUBLE,
                     neighbourRank, _bufferTag, _communicator, &request);

  if (result != MPI_SUCCESS) {
    std::cout << "rank " << _rank << " was not able to Irecv buffer from "
              << neighbourRank << ": " << std::endl;
    exit(EXIT_FAILURE);
  }

#if (MD_DEBUG == MD_YES)
  std::cout << "Rank " << _rank << " is receiving "
            << buffer->getLength() / 3 / MD_DIM << " molecules from "
            << neighbourRank << std::endl;
#endif
}
#endif
