// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder
#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/datastructures/FlexibleCellContainer.h"
#include "coupling/indexing/IndexingService.h"
#include <cstddef>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

namespace coupling {
namespace datastructures {
template <unsigned int dim> class FlexibleCellContainer;
}
} // namespace coupling

/**
 *	@brief provides access to coupling cells, which may belong to different indexing domains
 *
 * The FlexibleCellContainer is intended for use in cell domains which do not neatly correspond to a predefined indexing domain. For example, a CellContainer
 *cannot store cells in the macro2md overlap layers, and hence a FlexibleCellContainer is needed. The cells in a FlexibleCellContainer are assumed to be going
 *from one corner to another i.e. a cell is assumed to have a larger global scalar index than all of its preceding cells. Due to the "holes" expected due to the
 *non-contiguous nature of this container, index-based access is not possible or expected.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 */

template <unsigned int dim> class coupling::datastructures::FlexibleCellContainer {
public:
  FlexibleCellContainer() {}
  FlexibleCellContainer(std::vector<coupling::datastructures::CouplingCell<dim>*> couplingCells, std::vector<I01> idxs) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
    if (couplingCells.size() != idxs.size()) {
      std::cout << "ERROR size of index vector and coupling cell vector sent to FlexibleCellContainer constructor do not match";
      exit(EXIT_FAILURE);
    }
#endif
    _couplingCells.reserve(couplingCells.size());
    _idxs.reserve(idxs.size());
    for (std::size_t i = 0; i < couplingCells.size(); ++i) {
      _couplingCells.push_back(couplingCells[i]);
      _idxs.push_back(idxs[i]);
    }
  }

  template <class Container_T> FlexibleCellContainer(Container_T cells) {
    if constexpr (std::is_same_v<Container_T, FlexibleCellContainer>) {
      _idxs = cells._idxs;
      _couplingCells = cells._couplingCells;
    } else {
      auto numCells = cells.size();
      _idxs.reserve(numCells);
      _couplingCells.reserve(numCells);
      for (auto pair : cells)
        *this << pair;
    }
  }

  /** Adds a new coupling cell to the datastructure at the next index
   * @param cell a pointer to the cell to be inserted
   */
  void operator<<(std::pair<coupling::datastructures::CouplingCell<dim>*, I01> pair) {
    I01 idx;
    coupling::datastructures::CouplingCell<dim>* couplingCell;
    std::tie(couplingCell, idx) = pair;
    _couplingCells.push_back(couplingCell);
    _idxs.push_back(idx);
  }

  /**
   * Returns size of the underlying container.
   *
   * The number of indices stored should be equal to the number of cells stored, hence returning either is okay
   * @return the number of cells stored currently
   */
  unsigned int size() const { return _couplingCells.size(); }

  /**
   * @brief Provides iterator functionality (increment, access as <*cell, index> pair, equality)
   */
  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<dim>*>::const_iterator;
    using IndexIterator = std::vector<I01>::const_iterator;

    Iterator(CouplingCellIterator itCouplingCells, IndexIterator itIdxs) : _itCouplingCells(itCouplingCells), _itIdxs(itIdxs) {}

    /**
     * Iterator access, returning the data at the current iterator location
     *
     * @return a std::pair with the cell pointer and the index of the data that the iterator points to
     */
    const std::pair<coupling::datastructures::CouplingCell<dim>*, I01> operator*() const { return std::make_pair(*_itCouplingCells, *_itIdxs); }

    Iterator& operator++() {
      ++_itCouplingCells;
      ++_itIdxs;
      return *this;
    }

    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const Iterator& a, const Iterator& b) { return a._itCouplingCells == b._itCouplingCells && a._itIdxs == b._itIdxs; }

    friend bool operator!=(const Iterator& a, const Iterator& b) { return !(a == b); }

  private:
    /**Iterator to underlying cell* vector */
    CouplingCellIterator _itCouplingCells;

    /**Iterator to underlying index vector */
    IndexIterator _itIdxs;
  };
  /** Provides pointer to beginning of iterator of this container */
  Iterator begin() const { return Iterator(_couplingCells.begin(), _idxs.begin()); }

  /** Provides pointer to end of iterator of this container */
  Iterator end() const { return Iterator(_couplingCells.end(), _idxs.end()); }

private:
  /**Vector to store pointers to cells */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;

  /**Vector to store indices corresponding to cells in container
   * Since the cells can be part of any domain, they're stores as I01 (global noghost) type.
   */
  std::vector<I01> _idxs;
};