// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_LINKEDCELL_H_
#define _MOLECULARDYNAMICS_LINKEDCELL_H_

#include <list>
#include <iostream>
#include <cstdlib>
#include "simplemd/Molecule.h"

namespace simplemd { class LinkedCell; }

/** describes a linked cell.
 *  @author Philipp Neumann
 */
class simplemd::LinkedCell {
public:
  /** initialises linked cell list with numberMolecules empty entries */
  LinkedCell(const unsigned int numberMolecules = 0)
      : _molecules(numberMolecules) {}
  ~LinkedCell() { _molecules.clear(); }

  /** iterators to begin and end position */
  std::list<Molecule *>::iterator begin() { return _molecules.begin(); }
  std::list<Molecule *>::const_iterator constBegin() const {
    return _molecules.begin();
  }
  std::list<Molecule *>::iterator end() { return _molecules.end(); }
  std::list<Molecule *>::const_iterator constEnd() const {
    return _molecules.end();
  }

  std::list<Molecule *> &getList() { return _molecules; }
  const std::list<Molecule *> &getConstList() { return _molecules; }

  /** adds a molecule pointer */
  void addMolecule(Molecule *molecule) { _molecules.push_back(molecule); }

  /** deletes the molecule 'molecule' from the list, if it is contained */
  void deleteMolecule(Molecule *molecule) {
#if (MD_ERROR == MD_YES)
    if (molecule == NULL) {
      std::cout << "ERROR simplemd::LinkedCell::deleteMolecule: molecule==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    _molecules.remove(molecule);
  }

  /** initialises the molecule pointer at it with the value of molecule */
  void setMolecule(std::list<Molecule *>::iterator &it, Molecule *molecule) {
    *it = molecule;
  }

private:
  std::list<Molecule *> _molecules;
};

#endif // _MOLECULARDYNAMICS_LINKEDCELL_H_
