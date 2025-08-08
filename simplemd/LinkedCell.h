#ifndef _MOLECULARDYNAMICS_LINKEDCELL_H_
#define _MOLECULARDYNAMICS_LINKEDCELL_H_

#include "simplemd/Molecule.h"
#include <Kokkos_Core.hpp>
#include <iostream>

namespace simplemd {
class LinkedCell;
}

/**
 * A linked cell of a SimpleMD simulation for iterating the corresponding molecules.
 * NOTE: Never hold on to an instance of this class, as they are short lived contextual objects (host vs. device memory).
 * Always obtain them through simplemd::MoleculeContainer::operator[](int cellIndex).
 */
class simplemd::LinkedCell {
public:
  KOKKOS_FUNCTION LinkedCell() : _moleculeData(nullptr), _linkedCellNumMolecules(nullptr), _cellIndex(0), _isGhostCell(false) {}
  KOKKOS_FUNCTION LinkedCell(const Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace>* moleculeData,
                             const Kokkos::View<size_t*, Kokkos::LayoutRight, Kokkos::SharedSpace>* nMolecules, unsigned int cellIndex, bool isGhostCell)
      : _moleculeData(moleculeData), _linkedCellNumMolecules(nMolecules), _cellIndex(cellIndex), _isGhostCell(isGhostCell) {}

  class Iterator {
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Molecule;
    using pointer = value_type*;
    using reference = value_type&;

  public:
    KOKKOS_FUNCTION Iterator(pointer ptr) : _myPtr(ptr), _idx(0) {}

    KOKKOS_INLINE_FUNCTION reference operator*() const { return *_myPtr; }
    KOKKOS_INLINE_FUNCTION pointer operator->() { return _myPtr; }
    KOKKOS_INLINE_FUNCTION Iterator& operator++() {
      _myPtr++;
      _idx++;
      return *this;
    }
    KOKKOS_INLINE_FUNCTION Iterator operator++(int) {
      Iterator temp = *this;
      ++(*this);
      _idx++;
      return temp;
    }
    KOKKOS_INLINE_FUNCTION Iterator& operator--() {
      _myPtr--;
      _idx--;
      return *this;
    }
    KOKKOS_INLINE_FUNCTION Iterator operator--(int) {
      Iterator temp = *this;
      --(*this);
      _idx--;
      return temp;
    }

    KOKKOS_INLINE_FUNCTION friend bool operator==(const Iterator& a, const Iterator& b) { return a._myPtr == b._myPtr; }
    KOKKOS_INLINE_FUNCTION friend bool operator!=(const Iterator& a, const Iterator& b) { return a._myPtr != b._myPtr; }

    KOKKOS_INLINE_FUNCTION unsigned int getIndex() const { return _idx; }

  private:
    pointer _myPtr;
    unsigned int _idx;
  };

  KOKKOS_FUNCTION Iterator begin() const { return Iterator(getMolecule(0)); }
  KOKKOS_FUNCTION Iterator end() const { return Iterator(getMolecule(numMolecules())); }

  KOKKOS_FUNCTION void insert(Molecule& molecule) {
    *getMolecule(numMolecules()) = molecule;
    changeMoleculeCount(+1);
  }
  KOKKOS_FUNCTION void remove(int moleculeIdx) {
    *getMolecule(moleculeIdx) = *getMolecule(numMolecules() - 1);
    changeMoleculeCount(-1);
  }
  KOKKOS_FUNCTION void clear() { (*_linkedCellNumMolecules)(_cellIndex) = 0; }

  std::string to_string() const {
    std::stringstream to_ret;
    to_ret << " numMol: " << numMolecules() << std::endl;
    return to_ret.str();
  }

  KOKKOS_INLINE_FUNCTION unsigned int numMolecules() const { return (*_linkedCellNumMolecules)(_cellIndex); }

  KOKKOS_INLINE_FUNCTION bool isGhostCell() const { return _isGhostCell; }

  KOKKOS_INLINE_FUNCTION unsigned int getIndex() const { return _cellIndex; }

private:
  KOKKOS_INLINE_FUNCTION void changeMoleculeCount(int by) { (*_linkedCellNumMolecules)(_cellIndex) += by; }

  KOKKOS_INLINE_FUNCTION Molecule* getMolecule(unsigned int moleculeIndex) const { return &(*_moleculeData)(_cellIndex, moleculeIndex); }

  const Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace>* _moleculeData;
  const Kokkos::View<size_t*, Kokkos::LayoutRight, Kokkos::SharedSpace>* _linkedCellNumMolecules;
  const unsigned int _cellIndex;
  const bool _isGhostCell;
};

#endif // _MOLECULARDYNAMICS_LINKEDCELL_H_
