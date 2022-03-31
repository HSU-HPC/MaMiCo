// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MARDYNMOLECULEITERATOR_H_
#define MARDYNMOLECULEITERATOR_H_

#include "molecules/Molecule.h"

#include "coupling/interface/MoleculeIterator.h"
#include "coupling/interface/impl/MarDyn/MarDynMoleculeWrapper.h"
#include "coupling/interface/impl/MarDyn/MarDynCell.h"

/*	Implementation of the MoleculeIterator interface for Mardyn
 * 	@author Hanno Flohr, Philipp Neumann
 */
class MarDynMoleculeIterator
    : public coupling::interface::MoleculeIterator<MarDynCell, 3> {
public:
  /*	Constructor setting the iterator for the input MarDynCell,
 	 * 	the cutoff radius to the value of the input cell and
 	 * 	the MarDynMolecule '_mdMolecule' to NULL (prevents memory leaks in
 get-methods)
 	 */
  MarDynMoleculeIterator(MarDynCell &cell)
      : coupling::interface::MoleculeIterator<MarDynCell, 3>(cell),
        _cutoffRadius(cell.getCutoffRadius()) {}

  virtual ~MarDynMoleculeIterator() {}

  /* sets iterator to the first element of the cell */
  virtual void begin() {
    coupling::interface::MoleculeIterator<MarDynCell, 3>::_cell.begin();
  }

  /* returns true if the iterator should continue molecule traversal */
  virtual bool continueIteration() const {
    return coupling::interface::MoleculeIterator<MarDynCell, 3>::_cell
        .continueIteration();
  }

  /* sets iterator to next molecule of the cell */
  virtual void next() {
    coupling::interface::MoleculeIterator<MarDynCell, 3>::_cell.next();
  }

  /* returns a reference to the molecule that the iterator currently points to
   */
  virtual coupling::interface::Molecule<3> &get() {
    //set the mardyn molecule to the current iterator value
    _molecule.setMolecule(
        coupling::interface::MoleculeIterator<MarDynCell, 3>::_cell.get(),
        _cutoffRadius);
    return _molecule;
  }

  /* returns a const reference to the current molecule for pure reading purposes
   */
  virtual const coupling::interface::Molecule<3> &getConst() {
    //set the mardyn molecule to the current iterator value
    _molecule.setMolecule(
        coupling::interface::MoleculeIterator<MarDynCell, 3>::_cell.get(),
        _cutoffRadius);
    return _molecule;
  }

private:
  //the molecule wrapper interface returned in the get methods
  MarDynMoleculeWrapper _molecule;
  //the cutoff radius of the MarDyn cell
  const double _cutoffRadius;
};

#endif /* MARDYNMOLECULEITERATOR_H_ */
