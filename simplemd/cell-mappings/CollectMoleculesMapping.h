// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_COLLECTMOLECULESMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_COLLECTMOLECULESMAPPING_H_

#include "simplemd/LinkedCell.h"
#include "simplemd/services/MoleculeService.h"
#include <list>

namespace simplemd {
namespace cellmappings {
class CollectMoleculesMapping;
}
} // namespace simplemd

/** collects molecules from a linked cell, adds them to a buffer and removes
 * them from the linked cell list. This is used in the boundary treatment.
 *  @author Philipp Neumann
 */
class simplemd::cellmappings::CollectMoleculesMapping {
public:
  CollectMoleculesMapping(simplemd::services::MoleculeService &moleculeService) : _moleculeService(moleculeService) {}
  ~CollectMoleculesMapping() { _molecules.clear(); }

  void beginCellIteration() {}
  void endCellIteration() {}

  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    // loop over molecules
    for (std::list<Molecule *>::iterator it = cell.begin(); it != cell.end(); it++) {
      // push back molecule on buffer and remove it from simulation
      Molecule *myMolecule = (*it);
#if (MD_DEBUG == MD_YES)
      std::cout << "Delete molecule " << myMolecule->getConstPosition() << ", " << myMolecule->getConstVelocity() << " from MD" << std::endl;
#endif
      _molecules.push_back(*myMolecule);
      (*it) = NULL;

      // delete molecule from MoleculeService
      _moleculeService.deleteMolecule(*myMolecule);
    }
    cell.getList().clear();
  }
  void handleCellPair(LinkedCell &cell1, LinkedCell &cell2, const unsigned int &cellIndex1, const unsigned int &cellIndex2) {}

  void reset() { _molecules.clear(); }
  std::list<simplemd::Molecule> getCollectedMolecules() const { return _molecules; }

private:
  simplemd::services::MoleculeService &_moleculeService;
  std::list<simplemd::Molecule> _molecules;
};

#endif // _MOLECULARDYNAMICS_CELLMAPPINGS_COLLECTMOLECULESMAPPING_H_
