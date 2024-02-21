// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <vector>

#define DEBUG_JUNCTOR_INTERFACE

namespace coupling {
namespace filtering {
template <unsigned int dim, std::size_t inputc, std::size_t outputc> class JunctorInterface;
}
} // namespace coupling

/**
 * Junctors are a generalization of Filters, in the sense that they can harbour
 * multi in- and/or outputs. This allows for filters needing multiple sets of
 * input data (e.g. unfiltered/prefiltered) or comparing read-only-filters.
 * //TODO: implement FilterJunctorReadOnly? Junctors are stored in
 * FilterJunctions, which generalize FilterSequences in a similar way. You can
 * currently not add Junctors dynamically via FFF.
 *
 * Implemenents FI. The underlying FI has the junctor's main partition's (i.e
 * [0] of the corresponding std::array) input/output data. (cf. lines 35, 67.)
 * @author Felix Maurer
 */
template <unsigned int dim, std::size_t inputc, std::size_t outputc>
class coupling::filtering::JunctorInterface : public coupling::filtering::FilterInterface<dim> {
public:
  JunctorInterface(const std::array<std::vector<coupling::datastructures::CouplingCell<dim>*>, inputc> inputCellVectors,
                   const std::array<std::vector<coupling::datastructures::CouplingCell<dim>*>, outputc> outputCellVectors,
                   const std::array<bool, 7> filteredValues, const char* type)
      : // This assumes the array of cell vectors to be nonempty. Suboptimal.
        // NOTE: Avoid using FI's cell vectors. Use
        // _inputCellVectors/_outputCellVectors instead.
        coupling::filtering::FilterInterface<dim>(inputCellVectors[0], outputCellVectors[0], filteredValues, type), _inputCellVectors(inputCellVectors),
        _outputCellVectors(outputCellVectors) {
#ifdef DEBUG_JUNCTOR_INTERFACE
    // check input partition sizes
    for (unsigned int i = 0; i < inputc; i++)
      std::cout << "			FJ: Size of input cell vector at index " << i << ": " << inputCellVectors[i].size() << std::endl;
    for (unsigned int i = 0; i < outputc; i++)
      std::cout << "			FJ: Size of output cell vector at index " << i << ": " << outputCellVectors[i].size() << std::endl;

#endif
  }

  void updateCellData(std::vector<coupling::datastructures::CouplingCell<dim>*> new_inputCellVectors[inputc],
                      std::vector<coupling::datastructures::CouplingCell<dim>*> new_outputCellVectors[outputc],
                      std::vector<tarch::la::Vector<dim, unsigned int>>& new_cellIndices) {
    std::cout << "		JI: Updating cell data." << std::endl;
    _inputCellVectors = new_inputCellVectors;
    _outputCellVectors = new_outputCellVectors;

    // Assumes the input c-style vectors to be nonempty. May be problematic.
    coupling::filtering::FilterInterface<dim>::updateCellData(new_inputCellVectors[0], new_outputCellVectors[0]);
  }

protected:
  /**
   * Unlike regular filters, junctors allow for multiple input- and output-sets
   */
  std::array<std::vector<coupling::datastructures::CouplingCell<dim>*>, inputc> _inputCellVectors;
  std::array<std::vector<coupling::datastructures::CouplingCell<dim>*>, outputc> _outputCellVectors;
};
