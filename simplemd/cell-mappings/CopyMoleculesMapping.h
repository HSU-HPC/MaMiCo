// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_COPYMOLECULESMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_COPYMOLECULESMAPPING_H_

#include "simplemd/LinkedCell.h"
#include <list>

namespace simplemd {
namespace cellmappings {
class CopyMoleculesMapping;
}
} // namespace simplemd

/** creates a copy of all molecules and stores the molecules in a std-list.
 *  The molecules are sorted according to the linked cells that they belong to.
 *  After creating copies of the molecules, they are deleted from the
 * LinkedCell-structure.
 *  This mapping is used for reorganisation of the memory within the
 * MoleculeService.
 *
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::CopyMoleculesMapping {
public:
  CopyMoleculesMapping() {}
  ~CopyMoleculesMapping() {}

  void beginCellIteration() { _molecules.clear(); }
  void endCellIteration() {}
  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    // append molecules to list. Doing so, the molecules are sorted w.r.t. their
    // linked cell structure
    for (std::list<Molecule *>::iterator it = cell.begin(); it != cell.end(); it++) {
      _molecules.push_back(*(*it));
      // reset pointer to this molecule in the linked cell
      (*it) = NULL;
    }

    // delete molecules from the cell
    cell.getList().clear();
  }

  std::list<simplemd::Molecule> &getCopyOfMolecules() { return _molecules; }
  void removeCopy() { _molecules.clear(); }

private:
  std::list<simplemd::Molecule> _molecules;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_COPYMOLECULESMAPPING_H_
