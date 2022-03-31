// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _DUMMY_SOLVER_INTERFACE_SERVICE_H
#define _DUMMY_SOLVER_INTERFACE_SERVICE_H

#include "DummySolver.h"
#include "DummySolverInterface.h"
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include <cmath>

class DummySolverInterfaceService {
public:

  void init(tarch::la::Vector<3, unsigned int> numberProcesses,
            unsigned int rank, tarch::la::Vector<3, double> globalMDDomainSize,
            tarch::la::Vector<3, double> globalMDDomainOffset,
            tarch::la::Vector<3, double> macroscopicCellSize);

  void shutdown();

  static DummySolverInterfaceService &getInstance() {
    static DummySolverInterfaceService singleton;
    return singleton;
  }

  DummySolverInterface *getInterface() { return _dummySolverInterface; }

  /** add values to send buffer */
  bool addToSendBuffer(const double &density,
                       const tarch::la::Vector<3, double> &velocity,
                       const tarch::la::Vector<3, unsigned int> &index);

  /** receives values from receive buffer */
  bool getFromReceiveBuffer(
      double &density, tarch::la::Vector<3, double> &velocity,
      const tarch::la::Vector<3, unsigned int> &index) const;

  const unsigned int *getGlobalCellIndices4SendBuffer() const {
    return _globalIndices4SendBuffer;
  }
  unsigned int *getGlobalCellIndices4ReceiveBuffer() {
    return _globalIndices4ReceiveBuffer;
  }
  const std::vector<coupling::datastructures::MacroscopicCell<3> *> &
  getSendBuffer() const {
    return _sendBuffer;
  }
  const std::vector<coupling::datastructures::MacroscopicCell<3> *> &
  getReceiveBuffer() {
    return _receiveBuffer;
  }

  void resetSendBufferCounter() { _sendBufferCounter = 0; }

private:
  DummySolverInterfaceService()
      : _dummySolverInterface(NULL), _sendBufferSize(0), _receiveBufferSize(0),
        _globalIndices4SendBuffer(NULL), _globalIndices4ReceiveBuffer(NULL) {}
  ~DummySolverInterfaceService() {}

  /** loops over all macroscopic cells in the transfer region and allocates
    * buffers for all respective cells
    *  that are used in recv/send-operations.
    */
  void allocateBuffers();

  void deleteBuffersAndIndices();

  /** Create an object of the DummySolverInterface */
  DummySolverInterface *_dummySolverInterface;

  /** Returns true if the index is inside the transfer region */
  bool isInsideTransferRegion(
      const tarch::la::Vector<3, unsigned int> &index) const {
    bool isInside = true;
    tarch::la::Vector<3, unsigned int> _transferDomainOffsetIndex(2);
    tarch::la::Vector<3, unsigned int> _transferDomainNumCells(0);
    for (unsigned int d = 0; d < 3; d++) {
      _transferDomainNumCells[d] = _globalNumberMacroscopicCells[d] + 2;
    }
    for (unsigned int d = 0; d < 3; d++) {
      isInside = isInside && (index[d] >= _transferDomainOffsetIndex[d]) &&
                 (index[d] <
                  _transferDomainOffsetIndex[d] + _transferDomainNumCells[d]);
    }
    return isInside;
  }

  /** Returns global cell index (MaMiCo index) from the local index (dummy
   * solver index) */
  tarch::la::Vector<3, unsigned int> globalCellIndexfromLocalIndex(
      const tarch::la::Vector<3, unsigned int> _localIndex) const {
    tarch::la::Vector<3, unsigned int> _globalIndex(0);
    bool flag = isInsideTransferRegion(_localIndex);
    tarch::la::Vector<3, unsigned int> _transferDomainOffsetIndex(2);
    if (flag == true) {
      for (unsigned int d = 0; d < 3; d++) {
        _globalIndex[d] = _localIndex[d] - _transferDomainOffsetIndex[d];
      }
      return _globalIndex;
    } else {
      std::cout << "ERROR this cell does not lie in the transfer region"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  /** Returns Linearised cell index from global vector index*/
  unsigned int linearIndexFromVectorIndex(
      const tarch::la::Vector<3, unsigned int> globalVectorIndex) const {
    unsigned int linearIndex =
        globalVectorIndex[0] +
        globalVectorIndex[1] * (_globalNumberMacroscopicCells[0] + 2) +
        globalVectorIndex[2] * (_globalNumberMacroscopicCells[0] + 2) *
            (_globalNumberMacroscopicCells[1] + 2);
    return linearIndex;
  }

  /** buffer for send operations. */
  std::vector<coupling::datastructures::MacroscopicCell<3> *> _sendBuffer;
  /** buffer containing the corresponding cell identifiers for sending. */
  unsigned int *_globalIndices4SendBuffer;
  /** buffer for receive operations. */
  std::vector<coupling::datastructures::MacroscopicCell<3> *> _receiveBuffer;
  /** buffer containing the corresponding cell identifiers for receiving. */
  unsigned int *_globalIndices4ReceiveBuffer;

  /** send buffer counter */
  unsigned int _sendBufferCounter;
  /** size of the send buffer. Is determined in allocateBuffers() during
   * initialisation. */
  unsigned int _sendBufferSize;
  /** size of the receive buffer. Is determined in allocateBuffers() during
   * initialisation. */
  unsigned int _receiveBufferSize;

  /** size of MD domain + one macroscopic cell layer */
  tarch::la::Vector<3, double> _transferDomainSize;
  /** offset of MD domain - one macroscopic cell size */
  tarch::la::Vector<3, double> _transferDomainOffset;

  /** Macroscopic Cell Size as a vector of doubles */
  tarch::la::Vector<3, double> _globalMacroscopicCellSize;
  /** Number of macroscopic cells */
  tarch::la::Vector<3, unsigned int> _globalNumberMacroscopicCells;
};
#endif //_DUMMY_SOLVER_INTERFACE_SERVICE_H
