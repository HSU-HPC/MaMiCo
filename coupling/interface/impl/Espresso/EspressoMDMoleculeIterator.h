// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
#define _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_

#include "coupling/interface/MoleculeIterator.h"
#include "EspressoMDMolecule.h"
#include "particle_data.hpp"

namespace coupling {
namespace interface { class EspressoMDMoleculeIterator; }
}

/** molecule iterator.
 *  @author Rahul Arora
 */
class coupling::interface::EspressoMDMoleculeIterator
    : public coupling::interface::MoleculeIterator<ParticleList, 3> {
public:
  EspressoMDMoleculeIterator(ParticleList &cell)
      : coupling::interface::MoleculeIterator<ParticleList, 3>(cell),
        _buffer(NULL) {}

  virtual ~EspressoMDMoleculeIterator() {}

  /** sets the iterator to the first element */
  void begin() {
    _it = &_cell.part[0];
    _number = _cell.n;
    _counter = 0;
  }

  /** returns false, if the iterator reached the end of the molecule list */
  bool continueIteration() const { return (_counter < _number); }

  /** sets the iterator to the next molecule */
  void next() {
    // check if we need to update molecule information and reset flag in this
    // case
    _it++;
    _counter++;
  }

  /** returns a reference to the molecule that this iterator currently points to
   */
  coupling::interface::Molecule<3> &get() {
    _buffer.setMolecule(_it);
    return _buffer;
  }

  const coupling::interface::Molecule<3> &getConst() {
    _buffer.setMolecule(_it);
    return _buffer;
  }

private:
  int _number, _counter;
  Particle *_it;
  coupling::interface::EspressoMDMolecule _buffer;
};

#endif // _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
