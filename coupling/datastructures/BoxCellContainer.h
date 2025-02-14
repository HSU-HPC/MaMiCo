// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder
#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"

namespace coupling {
namespace datastructures {
class BoxCellContainer;
} // namespace datastructures
} // namespace coupling

/**
 * @brief provides access to the coupling cells. Base class for the class coupling::datastructures::LinkedCellContainer
 *
 * The BoxCellContainer is intended to be used for an arbitrary rectangular indexing domain defined by a lower and upper bound or lower bound and size.
 * Once created, it is assumed to hold all cells in that domain. BoxCellContainers do not capture relevant cells automatically.
 * It is the responsibility of the calling function to allocate the cells and populate the container completely before
 * using it. Due to this complete and contiguous nature of the datastructure, direct indexing is allowed, and the cell container is expected to start at
 * the relative location {0, 0, 0}. BoxCellContainers contain pointers to coupling cells, the cells are not owned by the container, hence
 * the calling function must delete the cells later.
 */
class coupling::datastructures::BoxCellContainer {

public:
  BoxCellContainer(I01 lowerBound, I01 upperBound) : BoxCellContainer({}, lowerBound, upperBound) {}
  
  BoxCellContainer(I01 lowerBound, tarch::la::Vector<3, int> shape) : BoxCellContainer({}, lowerBound, shape) {}

  BoxCellContainer(std::vector<coupling::datastructures::CouplingCell<3>*> couplingCells, I01 lowerBound, I01 upperBound)
    : BoxCellContainer(couplingCells, lowerBound, upperBound.get() - lowerBound.get()) {}

  BoxCellContainer(std::vector<coupling::datastructures::CouplingCell<3>*> couplingCells, I01 lowerBound, tarch::la::Vector<3, int> shape) {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    for (int i = 0; i < 3; i++)
    if (shape[i] < 0) {
        std::cout << "BoxCellContainer shape must not be negative, but dimension " << i << " was " << shape[i] << std::endl;
        std::exit(EXIT_FAILURE);       
    }
#endif
    _lowerBound = lowerBound;
    _shape = shape;
    _couplingCells.reserve(linearNumberCellsInDomain());
    for (auto cell : couplingCells) {
      *this << cell;
    }
  }

  /** Index based access, returns a pointer to the coupling cell
   * @param index index An index of type tarch::la::Vector
   * @return the pointer to the requested cell
   */
  coupling::datastructures::CouplingCell<3>* operator[](tarch::la::Vector<3, int> index) const {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    for (int i = 0; i < 3; i++)
    if (index[i] < 0 || index[i] >= numberCellsInDomain()[i]) {
      std::cout << "BoxCellContainer accessed out of bounds in dimension " << i;
      std::cout << ": " << index[i] << " not in 0-" << numberCellsInDomain()[i] << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (_couplingCells.size() < linearNumberCellsInDomain()) {
      std::cout << "BoxCellContainer accessed but not full " << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
    return _couplingCells[indexing::convertToScalar(_lowerBound + index)];
  }

  /** Adds a new coupling cell to the datastructure at the next index (will only work if the data structure is not yet full)
   * @param cell a pointer to the cell to be inserted
   */
  void operator<<(coupling::datastructures::CouplingCell<3>* couplingCell) {
#if (COUPLING_MD_ERROR == COUPLING_MD_YES)
    if (_couplingCells.size() >= linearNumberCellsInDomain()) {
      std::cout << "BoxCellContainer can only hold " << linearNumberCellsInDomain() << " coupling cells!"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
    _couplingCells.push_back(couplingCell);
  }

  /**
   * Returns size of the underlying container.
   *
   * This size is notably not the number of cells in the box domain. However, in practice the two numbers should be identical, as accessing contents of an
   * incomplete cotnainer is not allowed.
   * @return the number of cells stored currently
   */
  unsigned int size() const { return _couplingCells.size(); }

  /**
   * Returns shape of the underlying box domain.
   *
   * This size is notably not the number of cells currently held by the datastructure. However, in practice the two numbers should be identical, as accessing contents of an
   * incomplete cotnainer is not allowed.
   * @return the number of cells in the box domain along each dimension (x, y, z)
   */
  tarch::la::Vector<3, int> numberCellsInDomain() const {
    return _shape;
  }

  /**
   * Returns size of the underlying box domain.
   *
   * This size is notably not the number of cells currently held by the datastructure. However, in practice the two numbers should be identical, as accessing contents of an
   * incomplete cotnainer is not allowed.
   * @return the number of cells in the box domain
   */
  unsigned int linearNumberCellsInDomain() const {
    int count = 1;
    for (int i = 0; i < 3; i++)
        count *= _shape[i];
    return count;
  }

  // FIXME: Iterator does not do bounds check yet -> using wrong indices
  /**
   * @brief Provides iterator functionality (increment, access as <*cell, index> pair, equality)
   */
  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<3>*>::const_iterator;
    Iterator(CouplingCellIterator itCouplingCells, typename I01::IndexIterator itIdx) : _itCouplingCells(itCouplingCells), _itIdx(itIdx) {}
    /**
     * Iterator access, returning the data at the current iterator location
     *
     * @return a std::pair with the cell pointer and the index of the data that the iterator points to
     */
    const std::pair<coupling::datastructures::CouplingCell<3>*, I01> operator*() const { return std::make_pair(*_itCouplingCells, *_itIdx); }
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
    typename I01::IndexIterator _itIdx;
  };
  /** Provides pointer to beginning of iterator of this container */
  Iterator begin() const { return Iterator(_couplingCells.begin(), _lowerBound); }
  /** Provides pointer to end of iterator of this container */
  Iterator end() const { return Iterator(_couplingCells.end(), _upperBound()); }

protected:
  /**
   * Holds pointers to all coupling cells.
   */
  std::vector<coupling::datastructures::CouplingCell<3>*> _couplingCells;
  /**
   * Defines the lower corner of the box
   */
  I01 _lowerBound;
  /**
   * Defines the extent of the box along all dimensions
   */
  tarch::la::Vector<3, int> _shape;
  /**
   * Defines the upper corner of the box
   */
  I01 _upperBound() const {
    return I01(_lowerBound + _shape);
  }
  
};
