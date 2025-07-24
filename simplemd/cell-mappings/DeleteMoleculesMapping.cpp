// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/DeleteMoleculesMapping.h"

void simplemd::cellmappings::DeleteMoleculesMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  for (auto it = cell.begin(_moleculeService); it != cell.end(); ++it) {
    // get the molecule pointer and set it to NULL within the cell
    Molecule* myMolecule = (*it);

    // delete molecule from MoleculeService
    _moleculeService.deleteMolecule(*myMolecule);
  }
  cell.clear(_moleculeService);
}
