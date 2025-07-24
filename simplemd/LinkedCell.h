// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_LINKEDCELL_H_
#define _MOLECULARDYNAMICS_LINKEDCELL_H_

#include "simplemd/Molecule.h"
#include "simplemd/services/MoleculeService.h"
#include <cstdlib>
#include <iostream>
#include <list>

namespace simplemd {
class LinkedCell;
}

/** describes a linked cell.
 *  @author Philipp Neumann
 */
class simplemd::LinkedCell {
public:
  /** initialises linked cell list with numberMolecules empty entries */
  LinkedCell(const unsigned int index = -1) : _index(index), _moleculeCount(0) {}

  class Iterator {
  public:
    Iterator() : Iterator(0, 0) {} // For "uninitialized" stack variables in other classes
    Iterator(unsigned int cellIndex, unsigned int moleculeIndex, simplemd::services::MoleculeService* moleculeService = nullptr)
        : _moleculeService(moleculeService), _cellIndex(cellIndex), _moleculeIndex(moleculeIndex) {}
    Molecule* operator*() const {
      if (_moleculeService == nullptr)
        throw std::runtime_error("Cannot get molecule from simplemd::LinkedCell::Iterator (Pointer to simplemd::services::MoleculeService is a nullptr)");
      return _moleculeService->getCellMolecule(_cellIndex, _moleculeIndex);
    }
    Iterator operator++(int) {
      Iterator tmp = *this;
      _moleculeIndex++;
      return tmp;
    }
    Iterator operator--(int) {
      if (_moleculeIndex == 0)
        throw "Cannot decrement simplemd::LinkedCell::Iterator (_moleculeIndex cannot be negative)";
      Iterator tmp = *this;
      _moleculeIndex--;
      return tmp;
    }
    bool operator!=(const Iterator& other) const {
      if (_cellIndex != other._cellIndex)
        return true;
      return _moleculeIndex != other._moleculeIndex;
    }

  private:
    simplemd::services::MoleculeService* _moleculeService;
    unsigned int _cellIndex;
    unsigned int _moleculeIndex;
  };
  /** iterators to begin and end position */
  Iterator begin(simplemd::services::MoleculeService& moleculeService) const { return Iterator(_index, 0, &moleculeService); }
  Iterator end() const { return Iterator(_index, _moleculeCount); }

  void clear(simplemd::services::MoleculeService& moleculeService) {
    moleculeService.clearCellMolecules(_index);
    _moleculeCount = 0;
  }
  void addMolecule() { _moleculeCount++; }
  void deleteMolecule() {
    if (_moleculeCount == 0)
      throw std::runtime_error("Tried to call simplemd::LinkedCell::deleteMolecule() when _moleculeCount was already 0");
    _moleculeCount--;
  }
  unsigned int getMoleculeCount() const { return _moleculeCount; }
  unsigned int getIndex() const { return _index; }

private:
  unsigned int _index;
  unsigned int _moleculeCount;
};

#endif // _MOLECULARDYNAMICS_LINKEDCELL_H_
