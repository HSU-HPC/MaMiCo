// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/DeleteMoleculesMapping.h"

void simplemd::cellmappings::DeleteMoleculesMapping::
handleCell(LinkedCell& cell,const unsigned int &cellIndex) {
  // iterate over all molecules from this cell
  for (std::list<Molecule*>::iterator it = cell.begin(); it != cell.end(); it++){
    // get the molecule pointer and set it to NULL within the cell
    Molecule *myMolecule = (*it);
    (*it) = NULL;

    // delete molecule from MoleculeService
    _moleculeService.deleteMolecule(*myMolecule);
  }
  // clear cell list
  cell.getList().clear();
}
