// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/molecule-mappings/UpdateLinkedCellListsMapping.h"

void simplemd::moleculemappings::UpdateLinkedCellListsMapping::beginMoleculeIteration() {
  tarch::la::Vector<MD_DIM, unsigned int> bufferGlobal(0);
  tarch::la::Vector<MD_DIM, unsigned int> bufferLocal(0);
#if (MD_ERROR == MD_YES)
  _domainSize = _parallelTopologyService.getGlobalDomainSize();
#endif
  _domainOffset = _parallelTopologyService.getGlobalDomainOffset();
  _meshWidth = _parallelTopologyService.getMeshWidth();
  bufferGlobal = _parallelTopologyService.getGlobalIndexOfFirstCell();
  bufferLocal = _linkedCellService.getLocalIndexOfFirstCell();
  for (unsigned int d = 0; d < MD_DIM; d++) {
    _globalIndexOfFirstCell[d] = (int)bufferGlobal[d];
    _localIndexOfFirstCell[d] = (int)bufferLocal[d];
  }
}

void simplemd::moleculemappings::UpdateLinkedCellListsMapping::handleMolecule(Molecule& molecule) {

  const tarch::la::Vector<MD_DIM, double>& position(molecule.getConstPosition());
  for (unsigned int d = 0; d < MD_DIM; d++) {
#if (MD_ERROR == MD_YES)
    if ((position[d] < _domainOffset[d] - _meshWidth[d]) || (position[d] > _domainOffset[d] + _domainSize[d] + _meshWidth[d])) {
      std::cout << "ERROR "
                   "simplemd::moleculemappings::UpdateLinkedCellListsMapping::"
                   "handleMolecule: Position ";
      std::cout << d << " is out of range!" << std::endl;
      std::cout << "Position: " << position << std::endl;
      std::cout << molecule.getConstVelocity() << std::endl;
      std::cout << molecule.getConstForce() << std::endl;
      std::cout << molecule.getConstForceOld() << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
  }
  tarch::la::Vector<MD_DIM, unsigned int> cellVectorIndex(0);

  // determine current cell index (in serial, i.e. 1-D, form)
  for (unsigned int d = 0; d < MD_DIM; d++) {
    // find global cell index

    int index = (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d]));
    // shift into local cell index
    index += _localIndexOfFirstCell[d];
    index -= _globalIndexOfFirstCell[d];
#if (MD_ERROR == MD_YES)
    if (index < 0) {
      std::cout << "ERROR simplemd::moleculemappings::UpdateLinkedCellListsMapping: "
                   "index < 0: index=";
      std::cout << index << std::endl;
      std::cout << "Dimension : " << d << "," << _globalIndexOfFirstCell[d] << "," << _localIndexOfFirstCell[d] << std::endl;
      std::cout << (int)(floor((position[d] - _domainOffset[d]) / _meshWidth[d])) << std::endl;
      for (unsigned int e = 0; e < d; e++) {
        std::cout << cellVectorIndex[e] << std::endl;
      }
      std::cout << "Position: " << position << ", offset: " << _domainOffset << ", meshwidth: " << _meshWidth << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    cellVectorIndex[d] = (unsigned int)index;
  }

// if the molecule left the current cell, delete it from this cell and put it
// into the other cell
#if (MD_DEBUG == MD_YES)
  std::cout << "Put molecule " << molecule.getID() << ", pos. " << molecule.getConstPosition() << ", vel: " << molecule.getConstVelocity()
            << " , force: " << molecule.getConstForceOld() << " into cell " << std::endl;
  std::cout << cellVectorIndex << std::endl;
#endif
  _linkedCellService.addMoleculeToLinkedCell(molecule, cellVectorIndex);
}
