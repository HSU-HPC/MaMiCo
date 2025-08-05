#ifndef _MOLECULARDYNAMICS_LINKEDCELL_H_
#define _MOLECULARDYNAMICS_LINKEDCELL_H_

#include "simplemd/Molecule.h"
#include <Kokkos_Core.hpp>
#include <iostream>

namespace simplemd {
class LinkedCell;
}

class simplemd::LinkedCell {
public:
  KOKKOS_FUNCTION LinkedCell(Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> nMolecules,
                             Kokkos::View<Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> moleculeSlice)
      : _numMolecules(nMolecules), _moleculeData(moleculeSlice) {}

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

  KOKKOS_FUNCTION Iterator begin() const { return Iterator(&_moleculeData(0)); }
  KOKKOS_FUNCTION Iterator end() const { return Iterator(&_moleculeData(_numMolecules())); }

  KOKKOS_FUNCTION void insert(Molecule& molecule) {
    _moleculeData(_numMolecules()) = molecule;
    _numMolecules() += 1;
  }
  KOKKOS_FUNCTION void remove(int moleculeIdx) {
    _moleculeData(moleculeIdx) = _moleculeData(_numMolecules() - 1);
    _numMolecules() -= 1;
  }
  KOKKOS_FUNCTION void clear() { _numMolecules() = 0; }

  std::string to_string() const {
    std::stringstream to_ret;
    to_ret << " numMol: " << _numMolecules() << std::endl;
    return to_ret.str();
  }

  KOKKOS_INLINE_FUNCTION unsigned int numMolecules() const { return _numMolecules(); }

private:
  Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> _numMolecules;
  Kokkos::View<simplemd::Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> _moleculeData;
};

#endif // _MOLECULARDYNAMICS_LINKEDCELL_H_
