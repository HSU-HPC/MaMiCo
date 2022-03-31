// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "DummySolverInterfaceService.h"
#include <mpi.h>

void DummySolverInterfaceService::init(tarch::la::Vector<3, unsigned int> numberProcesses, unsigned int rank, tarch::la::Vector<3, double> globalMDDomainSize,
                                       tarch::la::Vector<3, double> globalMDDomainOffset, tarch::la::Vector<3, double> macroscopicCellSize) {
  _transferDomainSize = globalMDDomainSize + 2.0 * macroscopicCellSize;
  _transferDomainOffset = globalMDDomainOffset - macroscopicCellSize;
  _globalMacroscopicCellSize = macroscopicCellSize;
  for (unsigned int d = 0; d < 3; d++) {
    _globalNumberMacroscopicCells[d] = floor(globalMDDomainSize[d] / macroscopicCellSize[d] + 0.5);
    if (fabs(_globalNumberMacroscopicCells[d] * macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13) {
      std::cout << "ERROR MacroscopicSolverInterfaceService::init(): "
                   "globalNumberMacroscopicCells does not fit"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if (_dummySolverInterface != NULL) {
    delete _dummySolverInterface;
    _dummySolverInterface = NULL;
  }
  _dummySolverInterface = new DummySolverInterface(_globalNumberMacroscopicCells);
  if (_dummySolverInterface == NULL) {
    std::cout << "ERROR MacroscopicSolverInterfaceService::init(): Could not "
                 "initialise _transferStrategy!"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  allocateBuffers();
}

bool DummySolverInterfaceService::addToSendBuffer(const double &density, const tarch::la::Vector<3, double> &velocity,
                                                  const tarch::la::Vector<3, unsigned int> &index) {
  if (!isInsideTransferRegion(index)) {
    return false;
  }

  // check transfer strategy for send operations
  const tarch::la::Vector<3, unsigned int> globalVectorIndex = globalCellIndexfromLocalIndex(index);
  if (_dummySolverInterface->sendMacroscopicQuantityToMDSolver(globalVectorIndex)) {
    // if we have to send this cell, scale mass and momentum and write them into
    // cell buffer
    double volume = _globalMacroscopicCellSize[0] * _globalMacroscopicCellSize[1] * _globalMacroscopicCellSize[2];
    // volume = 1.0;
    const double mass = density * volume;
    const tarch::la::Vector<3, double> momentum = mass * velocity;
    _sendBuffer[_sendBufferCounter]->setMicroscopicMass(mass);
    _sendBuffer[_sendBufferCounter]->setMicroscopicMomentum(momentum);
    // store global index and increment send buffer-counter
    _globalIndices4SendBuffer[_sendBufferCounter] = linearIndexFromVectorIndex(globalVectorIndex);
    _sendBufferCounter++;
    return true;
  }
  return false;
}

bool DummySolverInterfaceService::getFromReceiveBuffer(double &density, tarch::la::Vector<3, double> &velocity,
                                                       const tarch::la::Vector<3, unsigned int> &index) const {
  if (!isInsideTransferRegion(index)) {
    return false;
  }

  // if we want to receive information of this cell...
  const tarch::la::Vector<3, unsigned int> globalVectorIndex = globalCellIndexfromLocalIndex(index);
  if (_dummySolverInterface->receiveMacroscopicQuantityFromMDSolver(globalVectorIndex)) {
    // ... search for cell in receive buffer...
    const unsigned int linearIndex = linearIndexFromVectorIndex(globalVectorIndex);
    double volume = _globalMacroscopicCellSize[0] * _globalMacroscopicCellSize[1] * _globalMacroscopicCellSize[2];
    // volume = 1.0;
    for (unsigned int i = 0; i < _receiveBufferSize; i++) {
      // ... and if we find the cell, write data to lbMass,lbMomentum and return
      if (_globalIndices4ReceiveBuffer[i] == linearIndex) {
        density = _receiveBuffer[i]->getMacroscopicMass() / volume;
        velocity = _receiveBuffer[i]->getMacroscopicMomentum() * (1.0 / (density * volume));
        return true;
      }
    }
    std::cout << "ERROR MacroscopicSolverInterfaceService::getFromReceiveBuffer(): "
                 "Cell information for global cell index "
              << globalVectorIndex << " not found in buffer!" << std::endl;
    exit(EXIT_FAILURE);
  }
  return false;
}

void DummySolverInterfaceService::allocateBuffers() {
  deleteBuffersAndIndices();
  // loop over global MD domain and find out if we need to receive/send a cell
  // store the cell indices for sending and receiving temporarily; we copy those
  // into the arrays afterwards, but first need to
  // determine the sizes of the arrays. For the send-buffer, this would
  // typically not be necessary, but we do this only once at
  // start-up, maybe it helps at debugging
  std::vector<unsigned int> receiveCellIndex;
  std::vector<unsigned int> sendCellIndex;
  tarch::la::Vector<3, unsigned int> loop(0);
  tarch::la::Vector<3, unsigned int> _transferDomainNumCells(0);
  for (unsigned int d = 0; d < 3; d++) {
    _transferDomainNumCells[d] = _globalNumberMacroscopicCells[d] + 2;
  }
  for (loop[2] = 0; loop[2] < _transferDomainNumCells[2]; loop[2]++) {
    for (loop[1] = 0; loop[1] < _transferDomainNumCells[1]; loop[1]++) {
      for (loop[0] = 0; loop[0] < _transferDomainNumCells[0]; loop[0]++) {
        const std::vector<unsigned int> ranks = _dummySolverInterface->getRanks(loop);
        const unsigned int sizeRanks = (unsigned int)ranks.size();
        bool cellOnThisProcess = false;
        for (unsigned int i = 0; i < sizeRanks; i++) {
          int thisRank;
          MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
          if (ranks[i] == thisRank) {
            cellOnThisProcess = true;
            break;
          }
        }

        // if the cell is on this process, find out if we need to receive/ send
        // information and increment buffers
        if (cellOnThisProcess) {
          if (_dummySolverInterface->receiveMacroscopicQuantityFromMDSolver(loop)) {
            receiveCellIndex.push_back(linearIndexFromVectorIndex(loop));
            _receiveBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
            if (_receiveBuffer[_receiveBufferSize] == NULL) {
              std::cout << "ERROR DummySolverInterfaceService::allocateBuffers(): "
                           "_receiveBuffer[i]==NULL!"
                        << std::endl;
              exit(EXIT_FAILURE);
            }
            _receiveBufferSize++;
          }

          if (_dummySolverInterface->sendMacroscopicQuantityToMDSolver(loop)) {
            sendCellIndex.push_back(linearIndexFromVectorIndex(loop));
            _sendBuffer.push_back(new coupling::datastructures::MacroscopicCell<3>());
            if (_sendBuffer[_sendBufferSize] == NULL) {
              std::cout << "ERROR DummySolverInterfaceService::allocateBuffers(): "
                           "_sendBuffer[i]==NULL!"
                        << std::endl;
              exit(EXIT_FAILURE);
            }
            _sendBufferSize++;
          }
        }
      }
    }
  } // loop over global domain

  // allocate indices
  if (_receiveBufferSize != 0) {
    // init global index buffer and copy cell indices into it
    _globalIndices4ReceiveBuffer = new unsigned int[_receiveBufferSize];
    if (_globalIndices4ReceiveBuffer == NULL) {
      std::cout << "ERROR DummySolverInterfaceService::allocateBuffers(): "
                   "_globalIndices4ReceiveBuffer==NULL"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < _receiveBufferSize; i++) {
      _globalIndices4ReceiveBuffer[i] = receiveCellIndex[i];
    }
    receiveCellIndex.clear();
  }

  if (_sendBufferSize != 0) {
    _globalIndices4SendBuffer = new unsigned int[_sendBufferSize];
    if (_globalIndices4SendBuffer == NULL) {
      std::cout << "ERROR DummySolverInterfaceService::allocateBuffers(): "
                   "_globalIndices4SendBuffer==NULL"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < _sendBufferSize; i++) {
      _globalIndices4SendBuffer[i] = sendCellIndex[i];
    }
    sendCellIndex.clear();
  }

  // reset counter for send buffer read-out
  _sendBufferCounter = 0;
}

void DummySolverInterfaceService::deleteBuffersAndIndices() {
  if (_globalIndices4SendBuffer != NULL) {
    delete[] _globalIndices4SendBuffer;
    _globalIndices4SendBuffer = NULL;
  }
  if (_globalIndices4ReceiveBuffer != NULL) {
    delete[] _globalIndices4ReceiveBuffer;
    _globalIndices4ReceiveBuffer = NULL;
  }
  for (unsigned int i = 0; i < _sendBufferSize; i++) {
    if (_sendBuffer[i] != NULL) {
      delete _sendBuffer[i];
      _sendBuffer[i] = NULL;
    }
  }
  _sendBuffer.clear();
  for (unsigned int i = 0; i < _receiveBufferSize; i++) {
    if (_receiveBuffer[i] != NULL) {
      delete _receiveBuffer[i];
      _receiveBuffer[i] = NULL;
    }
  }
  _receiveBuffer.clear();
  _sendBufferSize = 0;
  _receiveBufferSize = 0;
}

void DummySolverInterfaceService::shutdown() {
  if (_dummySolverInterface != NULL) {
    delete _dummySolverInterface;
    _dummySolverInterface = NULL;
  }
  deleteBuffersAndIndices();
}
