// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
#include <algorithm>
#include <vector>

namespace coupling {
namespace filtering {
template <class CellIndex_T, unsigned int dim> class Constant;
}
} // namespace coupling

/**
 * Filter applying a constant floating point value to all cells for all filtered
 * values. Optionally, one can specify to apply this filter only to certain
 * directions of multidimensional cell properties.
 *
 * @author Felix Maurer
 */

template <class CellIndex_T, unsigned int dim> class coupling::filtering::Constant : public coupling::filtering::FilterInterface<coupling::datastructures::CellContainer<CellIndex_T, dim>, dim> {
public:
  Constant(const coupling::datastructures::CellContainer<CellIndex_T, dim>& inputCellVector,
           const coupling::datastructures::CellContainer<CellIndex_T, dim>& outputCellVector, const std::array<bool, 7> filteredValues,
           const tarch::la::Vector<dim, bool> filteredDims, const double constant)
      : coupling::filtering::FilterInterface<coupling::datastructures::CellContainer<CellIndex_T, dim>, dim>(inputCellVector, outputCellVector, filteredValues, "Constant"), _constant(constant),
        _filteredDims(filteredDims) {}

  void operator()() {
    tarch::la::Vector<dim, double> vec_buf;
    for (unsigned int i = 0; i < this->_inputCells.size(); ++i) {
      // apply to scalars
      for (auto scalarProperty : this->_scalarAccessFunctionPairs) {
        (this->_outputCells[i]->*scalarProperty.set)(_constant);
      }

      // apply to vectors
      for (auto vectorProperty : this->_vectorAccessFunctionPairs) {
        // TODO: perhaps check if _filteredDims == true,..,true before this for
        // performance reasons?
        vec_buf = (this->_inputCells[i]->*vectorProperty.get)();

        for (unsigned int d = 0; d < dim; d++) {
          if (_filteredDims[d])
            vec_buf[d] = _constant;
        }

        (this->_outputCells[i]->*vectorProperty.set)(vec_buf);
      }
    }
  }

private:
  const double _constant;
  const tarch::la::Vector<dim, bool> _filteredDims;
};
