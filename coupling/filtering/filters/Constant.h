// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
#include <algorithm>
#include <vector>

namespace coupling {
namespace filtering {
template <unsigned int dim> class Constant;
}
} // namespace coupling

/**
 * Filter applying a constant floating point value to all cells for all filtered
 * values. Optionally, one can specify to apply this filter only to certain
 * directions of multidimensional cell properties.
 *
 * @author Felix Maurer
 */

template <unsigned int dim> class coupling::filtering::Constant : public coupling::filtering::FilterInterface<dim> {
public:
  Constant(const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& inputCellVector,
           const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& outputCellVector, const std::array<bool, 7> filteredValues,
           const tarch::la::Vector<dim, bool> filteredDims, const double constant)
      : coupling::filtering::FilterInterface<dim>(inputCellVector, outputCellVector, filteredValues, "Constant"), _constant(constant),
        _filteredDims(filteredDims) {}

  void operator()() {
    tarch::la::Vector<dim, double> vec_buf;
    for (auto outputCell : FilterInterface<dim>::_outputCells) {
      // apply to scalars
      for (auto scalarProperty : FilterInterface<dim>::_scalarAccessFunctionPairs) {
        (outputCell->*scalarProperty.set)(_constant);
      }

      // apply to vectors
      for (auto vectorProperty : FilterInterface<dim>::_vectorAccessFunctionPairs) {
        // TODO: perhaps check if _filteredDims == true,..,true before this for
        // performance reasons?
        vec_buf = (outputCell->*vectorProperty.get)();

        for (unsigned int d = 0; d < dim; d++) {
          if (_filteredDims[d])
            vec_buf[d] = _constant;
        }

        (outputCell->*vectorProperty.set)(vec_buf);
      }
    }
  }

private:
  const double _constant;
  const tarch::la::Vector<dim, bool> _filteredDims;
};
