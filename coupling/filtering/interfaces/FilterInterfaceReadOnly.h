// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "FilterInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class FilterInterfaceReadOnly;
}
} // namespace coupling

/*
 * Extension of FilterInterface.h for cases in which the filter itself does not
 * produce any output data. For such filters, you want to make use of
 * copyInputToOutput() (see below) in every filter step.
 *
 * @author Felix Maurer
 */

template <unsigned int dim> class coupling::filtering::FilterInterfaceReadOnly : public coupling::filtering::FilterInterface<dim> {
public:
  FilterInterfaceReadOnly(const std::vector<coupling::datastructures::MacroscopicCell<dim> *> &inputCellVector,
                          const std::vector<coupling::datastructures::MacroscopicCell<dim> *> &outputCellVector, const std::array<bool, 7> filteredValues,
                          const char *type)
      : coupling::filtering::FilterInterface<dim>(inputCellVector, outputCellVector, filteredValues, type) {}

protected:
  /*
   * Copies all filtered data from input to output. You always want to call this
   * as part of any implementation of coupling::FilterInterface<dim>::operator()
   * 	when implementing this interface, that is implementing a read-only
   * filter (e.g WriteToFile, Storuhal) If you would not do that, the successors
   * of the implementing filter in a sequence would get faulty input data.
   */
  void copyInputToOutput() {
    /*
     * In certain cases, e.g. read-only filters that operate on secondary cells
     * of an AsymmetricalFilterJunction, we don't want the filter to produce any
     * output. Then, and only then, the output cells vector will be empty.
     */
    if (coupling::filtering::FilterInterface<dim>::_outputCells.empty())
      return;

    for (unsigned int ci = 0; ci < coupling::filtering::FilterInterface<dim>::_outputCells.size(); ci++) {
      for (const auto scalarProperty : coupling::filtering::FilterInterface<dim>::_scalarAccessFunctionPairs) {
        (coupling::filtering::FilterInterface<dim>::_outputCells[ci]->*scalarProperty.set)(
            (coupling::filtering::FilterInterface<dim>::_inputCells[ci]->*scalarProperty.get)()); // call setter from output cell
                                                                                                  // using getter from input cell.
      }
      for (const auto vectorProperty : coupling::filtering::FilterInterface<dim>::_vectorAccessFunctionPairs) {
        (coupling::filtering::FilterInterface<dim>::_outputCells[ci]->*vectorProperty.set)(
            (coupling::filtering::FilterInterface<dim>::_inputCells[ci]->*vectorProperty.get)()); // call setter from output cell
                                                                                                  // using getter from input cell.
      }
    }
  }
};
