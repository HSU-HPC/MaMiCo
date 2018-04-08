// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEITERATOR_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEITERATOR_H_

#include "coupling/interface/Molecule.h"

namespace coupling {
  namespace interface {
    template<class LinkedCell,unsigned int dim>
    class MoleculeIterator;
  }
}


/** some iterator scheme for traversing the molecules within a linked cell.
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::interface::MoleculeIterator {
  protected:
    LinkedCell &_cell;
  public:
    MoleculeIterator(LinkedCell& cell): _cell(cell){}
    virtual ~MoleculeIterator(){}

    /** sets the iterator to the first element */
    virtual void begin() = 0;

    /** returns true, if the iterator should continue the molecule traversal */
    virtual bool continueIteration() const = 0;

    /** sets the iterator to the next molecule */
    virtual void next() = 0;

    /** returns a reference to the molecule that this iterator currently points to */
    virtual coupling::interface::Molecule<dim>& get() = 0;

    /** returns a const reference to the current molecule for pure reading purposes */
    virtual const coupling::interface::Molecule<dim>& getConst() = 0;
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEITERATOR_H_
