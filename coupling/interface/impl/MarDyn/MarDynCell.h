// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef MARDYNCELL_H_
#define MARDYNCELL_H_

#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"

/*	Helper class for the iteration over the elements of a Mardyn ParticleCell
 * 	@author Hanno Flohr
 */
class MarDynCell {
public:
  MarDynCell(ParticleCell *cell, double cutoffRadius) : _myCell(cell), _cutoffRadius(cutoffRadius) { _molecules = &(cell->getParticlePointers()); }

  MarDynCell() : _myCell(NULL), _molecules(NULL), _cutoffRadius(0.0) {}

  ~MarDynCell() {}

  // set iterator to first element
  void begin() { _moleculeIter = _molecules->begin(); }

  // increment the iterator
  void next() { _moleculeIter++; }

  // returns true if more molecules in the ParticleCell
  bool continueIteration() const { return (_moleculeIter != _molecules->end()); }

  // return the molecule the iterator is currently pointing to
  Molecule *get() { return *_moleculeIter; }

  // return the cutoff radius for this particle cell
  double getCutoffRadius() { return _cutoffRadius; }

  // return pointer to the particle cell
  ParticleCell *getParticleCell() { return _myCell; }

protected:
  ParticleCell *_myCell;
  std::vector<Molecule *> *_molecules;
  std::vector<Molecule *>::iterator _moleculeIter;
  double _cutoffRadius;
};

#endif /* MARDYNCELL_H_ */
