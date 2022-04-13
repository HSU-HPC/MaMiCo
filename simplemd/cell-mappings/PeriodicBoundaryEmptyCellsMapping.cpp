// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/PeriodicBoundaryEmptyCellsMapping.h"

simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::PeriodicBoundaryEmptyCellsMapping(
    simplemd::services::ParallelTopologyService& parallelTopologyService, simplemd::services::MoleculeService& moleculeService,
    simplemd::services::LinkedCellService& linkedCellService)
    : _parallelTopologyService(parallelTopologyService), _moleculeService(moleculeService), _linkedCellService(linkedCellService), _domainSize(0.0),
      _processCoordinates(0), _numberProcesses(0) {}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setDomainSize(const tarch::la::Vector<MD_DIM, double>& domainSize) { _domainSize = domainSize; }

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setProcessCoordinates(const tarch::la::Vector<MD_DIM, unsigned int>& processCoordinates) {
  _processCoordinates = processCoordinates;
}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::setNumberOfProcesses(const tarch::la::Vector<MD_DIM, unsigned int>& numberProcesses) {
  _numberProcesses = numberProcesses;
}

void simplemd::cellmappings::PeriodicBoundaryEmptyCellsMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  const tarch::la::Vector<MD_DIM, unsigned int> size(_linkedCellService.getLocalNumberOfCells() + 2u * _linkedCellService.getLocalIndexOfFirstCell());
  tarch::la::Vector<MD_DIM, unsigned int> outerCellCoords(0);
  tarch::la::Vector<MD_DIM, unsigned int> coords(0);
  simplemd::LinkedCell* innerCell = NULL;
  unsigned int helpIndex = cellIndex;
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
  for (std::list<Molecule*>::const_iterator it = cell.begin(); it != cell.end(); it++) {
    tarch::la::Vector<MD_DIM, double>& position = (*it)->getPosition();
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if ((outerCellCoords[d] == 0) && (_processCoordinates[d] == 0)) {
        position[d] += _domainSize[d];
      } else if ((outerCellCoords[d] == size[d] - 1) && (_processCoordinates[d] == _numberProcesses[d] - 1)) {
        position[d] -= _domainSize[d];
      }
    }
  }

  // if the molecules need to be sent, they are sent and deleted from the local molecule service
  if (_parallelTopologyService.reduceGhostCellViaBuffer(cell, cellIndex, _linkedCellService)) {
    for (std::list<Molecule*>::iterator it = cell.begin(); it != cell.end(); it++) {
      _moleculeService.deleteMolecule(*(*it));
    }
    // if the molecules need to be placed somewhere on this process, do so...
  } else {
    innerCell = &_linkedCellService.getLinkedCell(coords);
    // iterate over molecules and either send them to other process or put them locally in the right cell
    for (std::list<Molecule*>::iterator it = cell.begin(); it != cell.end(); it++) {
      Molecule myMolecule((*it)->getConstPosition(), (*it)->getConstVelocity());
      myMolecule.setForceOld((*it)->getConstForceOld());
      if ((*it)->isFixed())
        myMolecule.fix();
      _moleculeService.deleteMolecule(*(*it));
      Molecule* mPtr = _moleculeService.addMolecule(myMolecule);
      innerCell->addMolecule(mPtr);
    }
  }
  // in any case: clear list in this cell
  cell.getList().clear();
}
