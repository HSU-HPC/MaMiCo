// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder
#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"

namespace coupling {
namespace datastructures {
template <class CellIndexT, unsigned int dim> class CellContainer;
} // namespace datastructures
} // namespace coupling

/**
 * @brief provides access to the coupling cells. Base class for the class coupling::datastructures::LinkedCellContainer
 *
 * The CellContainer is intended to be used for a fixed indexing domain. Once created, it is assumed to hold all cells in that domain. For example, a
 * CellContainer created with template parameter I02 is expected to contain all cells that are accessible by iterating through I02. CellContainers do not
 * capture relevant cells automatically. It is the responsibility of the calling function to allocate the cells and populate the container completely before
 * using it. Due to this complete and contiguous nature of the datastructure, direct indexing is allowed, and the cell container is expected to start at
 * location 0 of the indexing domain it was initialized with. CellContainers contain pointers to coupling cells, the cells are not owned by the container, hence
 * the calling function must delete the cells later.
 * @tparam dim Number of dimensions; it can be 1, 2 or 3
 */
template <class CellIndexT, unsigned int dim> class coupling::datastructures::CellContainer {

public:
  CellContainer() { _couplingCells.reserve(CellIndexT::linearNumberCellsInDomain); }
  CellContainer(std::vector<coupling::datastructures::CouplingCell<dim>*> couplingCells) {
    _couplingCells.reserve(CellIndexT::linearNumberCellsInDomain);
    for (auto cell : couplingCells) {
      *this << cell;
    }
  }

  /** Index based access, returns a pointer to the coupling cell
   * @param index index An index of type CellIndexT i.e. from the initializing subdomain
   * @return the pointer to the requested cell
   */
  coupling::datastructures::CouplingCell<dim>* operator[](CellIndexT index) const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (_couplingCells.size() < CellIndexT::linearNumberCellsInDomain) {
      std::cout << "CellContainer<" << CellIndexT::TNAME << "," << dim << "> accessed but not full " << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
    return _couplingCells[indexing::convertToScalar(index)];
  }

  /** Adds a new coupling cell to the datastructure at the next index (will only work if the data structure is not yet full)
   * @param cell a pointer to the cell to be inserted
   */
  void operator<<(coupling::datastructures::CouplingCell<dim>* couplingCell) {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (_couplingCells.size() >= CellIndexT::linearNumberCellsInDomain) {
      std::cout << "CellContainer<" << CellIndexT::TNAME << "," << dim << "> can only hold " << CellIndexT::linearNumberCellsInDomain << " coupling cells!"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
    _couplingCells.push_back(couplingCell);
  }

  /**
   * Returns size of the underlying container.
   *
   * This size is notably not the number of cells in the index domain. However, in practice the two numbers should be identical, as accessing contents of an
   * incomplete cotnainer is not allowed.
   * @return the number of cells stored currently
   */
  unsigned int size() const { return _couplingCells.size(); }

  /**
   * @brief Provides iterator functionality (increment, access as <*cell, index> pair, equality)
   */
  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<dim>*>::const_iterator;

    Iterator(CouplingCellIterator itCouplingCells, typename CellIndexT::IndexIterator itIdx) : _itCouplingCells(itCouplingCells), _itIdx(itIdx) {}

    /**
     * Iterator access, returning the data at the current iterator location
     *
     * @return a std::pair with the cell pointer and the index of the data that the iterator points to
     */
    const std::pair<coupling::datastructures::CouplingCell<dim>*, CellIndexT> operator*() const { return std::make_pair(*_itCouplingCells, *_itIdx); }

    Iterator& operator++() {
      ++_itCouplingCells;
      ++_itIdx;
      return *this;
    }

    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const Iterator& a, const Iterator& b) { return a._itCouplingCells == b._itCouplingCells; }

    friend bool operator!=(const Iterator& a, const Iterator& b) { return !(a == b); }

  private:
    /**Iterator to underlying cell* vector */
    CouplingCellIterator _itCouplingCells;

    /**Iterator over the index subdomain
     *
     * As the container spans the whole index subdomain, the iterators already available for this index domain can be reused, and incremented in lockstep with
     * the internal iterator for the vector
     */
    typename CellIndexT::IndexIterator _itIdx;
  };

  /** Provides pointer to beginning of iterator of this container */
  Iterator begin() const { return Iterator(_couplingCells.begin(), CellIndexT::begin()); }

  /** Provides pointer to end of iterator of this container */
  Iterator end() const { return Iterator(_couplingCells.end(), CellIndexT::end()); }

protected:
  /**
   * Holds pointers to all coupling cells.
   * This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
};
