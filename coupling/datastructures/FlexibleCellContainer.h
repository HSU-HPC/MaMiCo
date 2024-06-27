// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#pragma once

#include <iterator>
#include <cstddef>
#include <vector>
#include <iostream>
#include <utility>
#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
#include "coupling/indexing/IndexingService.h"

namespace coupling {
namespace datastructures {
template <unsigned int dim> class FlexibleCellContainer;
}
}

template <unsigned int dim> class coupling::datastructures::FlexibleCellContainer {
public:

  void operator<<(std::pair<coupling::datastructures::CouplingCell<dim>*, I01*> pair) {
    I01* idx;
    coupling::datastructures::CouplingCell<dim>* couplingCell;
    std::tie(couplingCell, idx) = pair;
    _couplingCells.push_back(couplingCell);
    _idxs.push_back(idx);
  }

  int size() const {return _couplingCells.size();}

  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<dim>*>::iterator;
    using IndexIterator = std::vector<I01*>::iterator;

    Iterator(CouplingCellIterator itCouplingCells, IndexIterator itIdxs) : _itCouplingCells(itCouplingCells), _itIdxs(itIdxs) {}

    std::pair<coupling::datastructures::CouplingCell<dim>*, I01*> operator*() const { return std::make_pair(*_itCouplingCells, *_itIdxs); }

    Iterator& operator++() { return *this; }

    Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return *(a._itCouplingCells) == *(b._itCouplingCells) && *(a._itIdxs) == *(b._itIdxs); 
    }

    friend bool operator!=(const Iterator& a, const Iterator& b) { 
      return *(a._itCouplingCells) != *(b._itCouplingCells) || *(a._itIdxs) != *(b._itIdxs); 
    }

  private:
    CouplingCellIterator _itCouplingCells;
    IndexIterator _itIdxs;
  };

  Iterator begin() { return Iterator(_couplingCells.begin(), _idxs.begin()); }
  Iterator end()   { return Iterator(_couplingCells.end(), _idxs.end()); }
  

private:
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
  std::vector<I01*> _idxs;
};