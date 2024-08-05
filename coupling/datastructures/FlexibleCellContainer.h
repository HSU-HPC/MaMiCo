// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#pragma once

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/datastructures/CouplingCell.h"
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

  void operator<<(std::pair<coupling::datastructures::CouplingCell<dim>*, I01> pair) {
    I01 idx;
    coupling::datastructures::CouplingCell<dim>* couplingCell;
    std::tie(couplingCell, idx) = pair;
    _couplingCells.push_back(couplingCell);
    _idxs.push_back(idx);
  }

  int size() const { return _couplingCells.size(); }

  class Iterator {
  public:
    using CouplingCellIterator = typename std::vector<coupling::datastructures::CouplingCell<dim>*>::iterator;
    using IndexIterator = std::vector<I01>::iterator;

    Iterator(CouplingCellIterator itCouplingCells, IndexIterator itIdxs) : _itCouplingCells(itCouplingCells), _itIdxs(itIdxs) {}

    std::pair<coupling::datastructures::CouplingCell<dim>*, I01> operator*() const { return std::make_pair(*_itCouplingCells, _itIdxs); }

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

    friend bool operator==(const Iterator& a, const Iterator& b) { return *(a._itCouplingCells) == *(b._itCouplingCells) && a._itIdxs == b._itIdxs; }

    friend bool operator!=(const Iterator& a, const Iterator& b) { return *(a._itCouplingCells) != *(b._itCouplingCells) || a._itIdxs != b._itIdxs; }

  private:
    CouplingCellIterator _itCouplingCells;
    IndexIterator _itIdxs;
  };

  Iterator begin() { return Iterator(_couplingCells.begin(), _idxs.begin()); }
  Iterator end() { return Iterator(_couplingCells.end(), _idxs.end()); }

private:
  std::vector<coupling::datastructures::CouplingCell<dim>*> _couplingCells;
  std::vector<I01> _idxs;
};