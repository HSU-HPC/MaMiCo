// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
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
 *	@brief provides access to the coupling cells. Base class for the class
 *coupling::datastructures::LinkedCellContainer
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 */
template <class CellIndexT, unsigned int dim> class coupling::datastructures::CellContainer {

public:
  
  /** returns a pointer to the coupling cell without access to linked cells. */
  const coupling::datastructures::CouplingCell<dim>* operator[](CellIndexT index) const;

  /** adds a new coupling cell to the datastructure at the next index (will only work if the data structure is not yet full)*/
  void operator<<(coupling::datastructures::CouplingCell<dim>* couplingCell);

  int size() const {return _couplingCells.size();}

  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<dim>*>::iterator;

    Iterator(CouplingCellIterator itCouplingCells) : _itCouplingCells(itCouplingCells) {}

    std::pair<coupling::datastructures::CouplingCell<dim>*, I01*> operator*() const { 
      return std::make_pair(*_itCouplingCells, new CellIndexT(std::distance(_couplingCells.begin(), _itCouplingCells))); 
    }

    Iterator& operator++() { return *this; }

    Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return *(a._itCouplingCells) == *(b._itCouplingCells); 
    }

    friend bool operator!=(const Iterator& a, const Iterator& b) { 
      return *(a._itCouplingCells) != *(b._itCouplingCells); 
    }

  private:
    CouplingCellIterator _itCouplingCells;
  };

  Iterator begin() { return Iterator(_couplingCells.begin()); }
  Iterator end()   { return Iterator(_couplingCells.end()); }


protected:
  /** holds pointers to all coupling cells with linked cells, but without
   * access to linked cells. This is used for interfacing to send-recv
   * operations.
   */
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
};

#include "CellContainer.cpph"
