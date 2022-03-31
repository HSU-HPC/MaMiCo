// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
#define _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_

#include "simplemd/Molecule.h"
#include "simplemd/LinkedCell.h"
#include "coupling/interface/MoleculeIterator.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDMolecule.h"

namespace coupling {
namespace interface {

/** iterates over molecules in a SimpleMD linked cell.
 *  @author Philipp Neumann
 */
class SimpleMDMoleculeIterator
    : public MoleculeIterator<simplemd::LinkedCell, MD_DIM> {
public:
  SimpleMDMoleculeIterator(simplemd::LinkedCell &cell)
      : coupling::interface::MoleculeIterator<simplemd::LinkedCell, MD_DIM>(
            cell),
        _buffer(NULL) {}
  virtual ~SimpleMDMoleculeIterator() {}

  /** sets the iterator to the first element */
  void begin() { _it = _cell.begin(); }

  /** returns false, if the iterator reached the end of the molecule list */
  bool continueIteration() const { return (_it != _cell.end()); }

  /** sets the iterator to the next molecule */
  void next() {
    // check if we need to update molecule information and reset flag in this
    // case
    _it++;
  }

  /** returns a reference to the molecule that this iterator currently points to
   */
  coupling::interface::Molecule<MD_DIM> &get() {
    _buffer.setMolecule(*_it);
    return _buffer;
  }

  const coupling::interface::Molecule<MD_DIM> &getConst() {
    _buffer.setMolecule(*_it);
    return _buffer;
  }

private:

  std::list<simplemd::Molecule *>::iterator _it;
  coupling::interface::SimpleMDMolecule _buffer;
};

}
}

#endif // _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
