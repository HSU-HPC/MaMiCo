// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterfaceReadOnly.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class Copy;
}
} // namespace coupling

/**
 * Simplest practical implementation of FilterInterfaceReadOnly.
 * Used primarily as a "filler" filter when adding filters at runtime. Check
 * FilterSequence.h for details.
 * @author Felix Maurer
 */

template <unsigned int dim> class coupling::filtering::Copy : public coupling::filtering::FilterInterfaceReadOnly<dim> {
public:
  Copy(const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& inputCellVector,
       const std::vector<coupling::datastructures::MacroscopicCell<dim>*>& outputCellVector, const std::array<bool, 7> filteredValues)
      : coupling::filtering::FilterInterfaceReadOnly<dim>(inputCellVector, outputCellVector, filteredValues, "COPY") {}

  void operator()() {
    coupling::filtering::FilterInterfaceReadOnly<dim>::copyInputToOutput();
    // coupling::filtering::FilterInterface<dim>::DEBUG_PRINT_CELL_VELOCITY("COPY");
  }
};
