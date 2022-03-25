// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/PeriodicAndParallelBoundaryFillCellsMapping.h"

void simplemd::cellmappings::PeriodicAndParallelBoundaryFillCellsMapping::handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
  // size of the local domain
  const tarch::la::Vector<MD_DIM, unsigned int> size(_linkedCellService.getLocalNumberOfCells() + 2u * _linkedCellService.getLocalIndexOfFirstCell());
  // first: send molecules to neighboring ghost cells, if needed.
  std::vector<tarch::la::Vector<MD_DIM, unsigned int>> localIndex = _parallelTopologyService.broadcastInnerCellViaBuffer(cell, cellIndex, _linkedCellService);

  // now: run over the local periodic boundaries and update those
  const unsigned int localIndexSize = (unsigned int)localIndex.size();
  for (unsigned int i = 0; i < localIndexSize; i++) {
    // compute shift for periodic condition
    tarch::la::Vector<MD_DIM, double> shift(0.0);
    tarch::la::Vector<MD_DIM, double> positionBuffer(0.0);
    for (unsigned int d = 0; d < MD_DIM; d++) {
      if (localIndex[i][d] == 0) {
        shift[d] = -_domainSize[d];
      } else if (localIndex[i][d] == size[d] - 1) {
        shift[d] = _domainSize[d];
      }
    }
    // iterate over cell's molecules and add molecules to local ghost cell
    for (std::list<Molecule *>::const_iterator it = cell.begin(); it != cell.end(); it++) {
      positionBuffer = (*it)->getConstPosition();
      positionBuffer += shift;
      Molecule myMolecule(positionBuffer, (*it)->getConstVelocity());
      myMolecule.setForceOld((*it)->getConstForceOld());
      if ((*it)->isFixed())
        myMolecule.fix();
      Molecule *myPtr = _moleculeService.addMolecule(myMolecule);
      _linkedCellService.addMoleculeToLinkedCell(*myPtr, localIndex[i]);
    }
  }
}
