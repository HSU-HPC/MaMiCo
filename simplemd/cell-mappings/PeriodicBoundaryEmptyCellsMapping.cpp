// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.h"

simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::PeriodicBoundaryEmptyCellsMapping(
    simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::MoleculeContainer& moleculeContainer)
    : _parallelTopologyService(parallelTopologyService), _moleculeContainer(moleculeContainer), _domainSize(0.0), _processCoordinates(0), _numberProcesses(0) {}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setDomainSize(const tarch::la::Vector<MD_DIM, double>& domainSize) { _domainSize = domainSize; }

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setProcessCoordinates(const tarch::la::Vector<MD_DIM, unsigned int>& processCoordinates) {
  _processCoordinates = processCoordinates;
}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setNumberOfProcesses(const tarch::la::Vector<MD_DIM, unsigned int>& numberProcesses) {
  _numberProcesses = numberProcesses;
}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::handleCell(LinkedCell& cell) {
  const tarch::la::Vector<MD_DIM, unsigned int> size(_moleculeContainer.getLocalNumberOfCells() + 2u * _moleculeContainer.getLocalIndexOfFirstCell());
  tarch::la::Vector<MD_DIM, unsigned int> outerCellCoords(0);
  tarch::la::Vector<MD_DIM, unsigned int> coords(0);
  size_t helpIndex = cell.getIndex();
#if (MD_DIM > 2)
  coords[2] = helpIndex / (size[1] * size[0]);
  helpIndex = helpIndex - coords[2] * (size[1] * size[0]);
#endif
#if (MD_DIM > 1)
  coords[1] = helpIndex / size[0];
  helpIndex = helpIndex - coords[1] * size[0];
#endif
  coords[0] = helpIndex;
  // now, coords and outerCell are identical
  outerCellCoords = coords;

  // after this, coords contains the local coordinates of the inner cell on the respective process
  for (unsigned int d = 0; d < MD_DIM; d++) {
    if ((coords[d] == 0) && (_processCoordinates[d] == 0)) {
      coords[d] = size[d] - 2;
    } else if ((coords[d] == size[d] - 1) && (_processCoordinates[d] == _numberProcesses[d] - 1)) {
      coords[d] = 1;
    }
  }

  // now: iterate over all molecules within this cell and change position
  for (auto it = cell.begin(); it != cell.end(); it++) {
    tarch::la::Vector<MD_DIM, double>& position = it->getPosition();
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if ((outerCellCoords[d] == 0) && (_processCoordinates[d] == 0)) {
        position[d] += _domainSize[d];
      } else if ((outerCellCoords[d] == size[d] - 1) && (_processCoordinates[d] == _numberProcesses[d] - 1)) {
        position[d] -= _domainSize[d];
      }
    }
  }

  // if the molecules need to be sent, they are sent and deleted from the local molecule service
  if (!_parallelTopologyService.reduceGhostCellViaBuffer(cell, cell.getIndex(), _moleculeContainer)) {
    // if the molecules need to be placed somewhere on this process, do so...
    _moleculeContainer.sort(_moleculeContainer[outerCellCoords].getIndex());
  }
  cell.clear();
}
