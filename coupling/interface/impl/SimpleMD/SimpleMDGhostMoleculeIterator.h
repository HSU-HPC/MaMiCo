// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDGHOSTMOLECULEITERATOR_CPP_
#define _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDGHOSTMOLECULEITERATOR_CPP_

#include "simplemd/Molecule.h"
#include "simplemd/LinkedCell.h"
#include "coupling/interface/MoleculeIterator.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDMolecule.h"

namespace coupling {
namespace interface {

/** iterates over ghost molecules in a SimpleMD linked cell.
 *  @author Helene Wittenberg
 */
class SimpleMDGhostMoleculeIterator: public MoleculeIterator<simplemd::LinkedCell,MD_DIM> {
  public:
    SimpleMDGhostMoleculeIterator(simplemd::LinkedCell& cell):
    coupling::interface::MoleculeIterator<simplemd::LinkedCell,MD_DIM>(cell),
    _buffer(NULL){}
    virtual ~SimpleMDGhostMoleculeIterator(){}

    /** sets the iterator to the first element */
    void begin()override{
      _it = _cell.beginGhost();
    }

    /** sets the iterator to the first element */
    void begin(const unsigned int dim){
      _it = _cell.beginGhost(dim);
    }

    /** returns false, if the iterator reached the end of the molecule list */
    bool continueIteration() const override{
      return (_it != _cell.endGhost());
    }

    /** returns false, if the iterator reached the end of the molecule list */
    bool continueIteration(const unsigned int dim) const {
      return (_it != _cell.endGhost(dim));
    }

    /** sets the iterator to the next molecule */
    void next()override{
      _it++;
    }

    /** returns a reference to the molecule that this iterator currently points to */
    coupling::interface::Molecule<MD_DIM>& get(){
      _buffer.setMolecule(*_it);
      return _buffer;
    }

    const coupling::interface::Molecule<MD_DIM>& getConst(){
      _buffer.setMolecule(*_it);
      return _buffer;
    }

  private:
    std::list<simplemd::Molecule*>::iterator _it;
    coupling::interface::SimpleMDMolecule _buffer;
};
}
}

#endif // _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDGHOSTMOLECULEITERATOR_CPP_
