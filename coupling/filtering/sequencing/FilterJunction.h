// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_JUNCTION

#include "coupling/filtering/sequencing/FilterSequence.h"

// INCLUDE ALL JUNCTOR HEADERS HERE
#include "coupling/filtering/filters/NLM.h"
#include "coupling/filtering/interfaces/JunctorInterface.h"

/*
 * Generalizes the concept of FilterSequences: A FilterJunction can have more
 * than a single input. This allows for seemless integration of FilterJunctions
 * into a net of FilterSequences.
 *
 *
 * @todo Support multiple outputs. ("X-Junctions")
 * @todo Support dynamically linked filters.
 * @author Felix Maurer
 */

namespace coupling {
namespace filtering {
template <class Container_T, unsigned int dim, std::size_t inputc> class FilterJunction;
}
} // namespace coupling

template <class Container_T, unsigned int dim, std::size_t inputc> class coupling::filtering::FilterJunction : public coupling::filtering::FilterSequence<Container_T, dim> {
public:
  FilterJunction(const char* name,
                 const std::vector<Container_T> inputCellVector, // inputc input cell containers
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                 MPI_Comm comm,
#endif
                 std::array<bool, 7> filteredValues)
      : coupling::filtering::FilterSequence<Container_T, dim>(name, inputCellVector[0],
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                 comm,
#endif
                                                 filteredValues) {
    if (inputc == 0)
      throw std::runtime_error("ERROR: Creating FilterJunction with inputc = 0.");

    if (inputc != inputCellVector.size())
      throw std::runtime_error("ERROR: Creating FilterJunction, wrong number of inputs");

#ifdef DEBUG_FILTER_JUNCTION
    std::cout << PRINT_PREFIX() << "This is a FilterJunction. Number of inputs:" << inputc << std::endl;
#endif

    for (unsigned int p = 0; p < inputc; p++) 
      _inputCellVector_parted[p] = inputCellVector[p];

    coupling::filtering::FilterSequence<Container_T, dim>::_isModifiable = false; // Dynamic filters are not yet supported.
  }

  /*
   * This member function allows appendance and insertion of filters defined by
   * two processing functions to a modifiable sequence at runtime. Index -1
   * implies appending.
   */
  void addFilter(
      const std::function<std::vector<double>(std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* applyScalar,
      const std::function<std::vector<std::array<double, dim>>(std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* applyVector,
      int filterIndex = -1) override {
// Do nothing, not yet supported.
#ifdef DEBUG_FILTER_JUNCTION
    std::cout << PRINT_PREFIX()
              << "This is a FilterJunction. addFilter(...) is not supported "
                 "and has no effect."
              << std::endl;
#endif
  }

  /*
   * This function is very similar to the interface's. Check
   * coupling::FilterSequence for more details.
   */
  int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode) override;

  /*
   * The first partition of _cellVector1/2 is the main partition. A junction's
   * default output is always its main partition.
   */
  const Container_T& getOutputCellVector(unsigned int outputIndex = 0) const override {
    if (outputIndex >= inputc) {
      std::cout << PRINT_PREFIX() << "ERROR: getOutputCellVector: Requested output index(" << outputIndex << ") too high. (partitions: )" << inputc
                << std::endl;
      exit(EXIT_FAILURE);
    }

    if (coupling::filtering::FilterSequence<Container_T, dim>::_filters.empty())
      std::cout << PRINT_PREFIX() << "Warning: Accessing cell vectors while _filters is empty." << std::endl;
    if (coupling::filtering::FilterSequence<Container_T, dim>::_filters.size() % 2 == 0)
      return this->_cellVector1;
    else
      return this->_cellVector2;
  }

  void printFilters() override {
    std::cout << "Junctors in junction " << coupling::filtering::FilterSequence<Container_T, dim>::_name << ": ";
    for (auto f : coupling::filtering::FilterSequence<Container_T, dim>::_filters)
      std::cout << f->getType() << " ";
    std::cout << std::endl;
  }

  std::string PRINT_PREFIX() const override {
    return std::string("	FJ(").std::string::append(coupling::filtering::FilterSequence<Container_T, dim>::_name).std::string::append("): ");
  }

private:
  // These must be parted for junction output and junctor in/output.
  std::array<Container_T, inputc> _inputCellVector_parted;
};

// inlcude implementation
#include "coupling/filtering/sequencing/FilterJunction.cpph"
