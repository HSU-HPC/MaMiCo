// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef LMP_MAMICO_CELL_H_
#define LMP_MAMICO_CELL_H_

#include <vector>

namespace LAMMPS_NS {
// implements the storage of molecule information in form of 2D pointer
// structures
class MoleculeInformation {
public:
  double **_x;
  double **_v;
  double **_f;
};

// implements a cell structure which embed atom identifiers
class MamicoCell {
public:
  MamicoCell() : _it(_atoms.begin()), _isGhostCell(false) {}
  ~MamicoCell() {}

  // deletes list of atom identifiers
  void clear() { _atoms.clear(); }
  // adds an atom id
  void addAtom(int id) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    std::cout << "MamicoCell::addAtom(): Add atom " << id << std::endl;
#endif
    _atoms.push_back(id);
  }

  // delete an atom id
  void deleteAtom(int id) {
    const unsigned int size = _atoms.size(); // number of atoms in this cell
    unsigned int index = size; // index which contains the atom to be deleted
    for (unsigned int i = 0; i < size; i++) {
      if (_atoms[i] == id)
        index = i;
    }
    // if atom is contained, delete it from vector; therefore, shift entries
    // towards front of vector and remove last vector entry
    if (index != size) {
      for (unsigned int i = index; i < size - 1; i++) {
        _atoms[i] = _atoms[i + 1];
      }
      _atoms.pop_back();
    }
  }

  // sets the iterator to first element registered atom ids
  void begin() { _it = _atoms.begin(); }
  // increments iterator
  void next() { _it++; }
  // returns true, if there are more atoms in this cell
  bool continueIteration() const { return (_it != _atoms.end()); }
  // returns the id of the atom at the current iterator position
  int get() { return *_it; }

  bool isGhostCell() const { return _isGhostCell; }
  void setGhostCell(bool isGhostCell) { _isGhostCell = isGhostCell; }

private:
  std::vector<int> _atoms;        // stores the atom ids in a vector
  std::vector<int>::iterator _it; // iterates over the vector
  bool _isGhostCell; // true, if this cell is flagged as a ghost cell
};
} // namespace LAMMPS_NS

#endif // LMP_MAMICO_CELL_H
