// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_JUNCTION_ASYM

#include "coupling/filtering/sequencing/FilterSequence.h"

// INCLUDE ALL RELEVANT JUNCTOR HEADERS HERE
#include "coupling/filtering/filters/WriteToFileJunctor.h"

/*
 * Filtering junction which does not use multiple identical input domains, but
 * rather uses one primary input partition to filter on while having access to a
 * secondary input data set. As no filtering is conducted on that secondary
 * input cell data set, it may be of any shape or size. The primary input
 * partition has to follow the restrictions to input data sets of
 * FilterSequences and FilterJunctions.
 *
 * @author Felix Maurer
 */

namespace coupling {
namespace filtering {
template <unsigned int dim> class AsymmetricalFilterJunction;
}
} // namespace coupling

template <unsigned int dim>
class coupling::filtering::AsymmetricalFilterJunction
    : public coupling::filtering::FilterSequence<dim> {
public:
  AsymmetricalFilterJunction(
      const char *name,
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
          primaryInputCellVector, // primary input of sequence.
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
          secondaryInputCellVector, // additional data, presented as macro cells
                                    // as well
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
      MPI_Comm comm,
#endif
      std::array<bool, 7> filteredValues)
      : coupling::filtering::FilterSequence<dim>(name, primaryInputCellVector,
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                 comm,
#endif
                                                 filteredValues),
        _inputCellVector_secondary(secondaryInputCellVector) {
#ifdef DEBUG_FILTER_JUNCTION_ASYM
    std::cout << PRINT_PREFIX() << "Begin initialization." << std::endl;
#endif

    // allocate and init secondary cell vectors
    for (auto cell : _inputCellVector_secondary) {
      _cellVector1_secondary.push_back(
          new coupling::datastructures::MacroscopicCell<dim>(*cell));
      _cellVector2_secondary.push_back(
          new coupling::datastructures::MacroscopicCell<dim>(*cell));
    }
#ifdef DEBUG_FILTER_JUNCTION_ASYM
    std::cout << PRINT_PREFIX() << "Initialized secondary cell vectors."
              << std::endl;
    std::cout << PRINT_PREFIX()
              << "First element of _cellVector1_secondary after init: "
              << _cellVector1_secondary[0] << std::endl;
    std::cout << PRINT_PREFIX()
              << "First element of _cellVector2_secondary after init: "
              << _cellVector2_secondary[0] << std::endl;
#endif

    coupling::filtering::FilterSequence<dim>::_isModifiable =
        false; // Dynamic filters are not yet supported. TODO
  }

  ~AsymmetricalFilterJunction() {
    for (auto secondarycell : _cellVector1_secondary)
      delete secondarycell;
    for (auto secondarycell : _cellVector2_secondary)
      delete secondarycell;
  }

  /*
   * This function is very similar to the interface's. Check
   * coupling::FilterSequence for more details.
   */
  int loadFiltersFromXML(tinyxml2::XMLElement *sequenceNode) override;

  void printFilters() override {
    std::cout << "Junctors in asymmetrical junction "
              << coupling::filtering::FilterSequence<dim>::_name << ": ";
    for (auto f : coupling::filtering::FilterSequence<dim>::_filters)
      std::cout << f->getType() << " ";
    std::cout << std::endl;
  }

  std::string PRINT_PREFIX() const override {
    return std::string("	AFJ(")
        .std::string::append(coupling::filtering::FilterSequence<dim>::_name)
        .std::string::append("): ");
  }

private:
  std::vector<coupling::datastructures::MacroscopicCell<dim> *>
      _inputCellVector_secondary;

  std::vector<coupling::datastructures::MacroscopicCell<dim> *>
      _cellVector1_secondary;
  std::vector<coupling::datastructures::MacroscopicCell<dim> *>
      _cellVector2_secondary;
};

// inlcude implementation
#include "coupling/filtering/sequencing/AsymmetricalFilterJunction.cpph"
