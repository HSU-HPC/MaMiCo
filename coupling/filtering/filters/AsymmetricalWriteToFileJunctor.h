// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_WTF_JUNCTION

#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/interfaces/AsymmetricalJunctorInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class AsymmetricalWriteToFileJunctor;
}
} // namespace coupling

/**
 * Takes the inner and outer MD domains and combines them to
 * one vector representing the entire domain. For reasons of
 * simplicity, first a dummy filter is used to ensure proper
 * junctor output. Next the real write to file is This vector is
 * then written to a single file.
 * @author Sebastian Lerdo
 */
template <unsigned int dim> class coupling::filtering::AsymmetricalWriteToFileJunctor : public coupling::filtering::AsymmetricalJunctorInterface<dim> {
public:
  AsymmetricalWriteToFileJunctor(
      // first cell data set
      const std::vector<coupling::datastructures::MacroscopicCell<dim>*> inputCellVector1,
      const std::vector<coupling::datastructures::MacroscopicCell<dim>*> outputCellVector1,
      // second cell data set
      const std::vector<coupling::datastructures::MacroscopicCell<dim>*> inputCellVector2,
      // no secondary output cells
      //"global" parameters for both WriteToFile instances
      const std::array<bool, 7> filteredValues,
      // WriteToFile-specific parameters. [0] is for the dummy WriteToFile instance, while [1] is for the desired WriteToFile instance
      std::array<std::string, 2> location, std::array<bool, 2> overwrite = {false}, std::array<int, 2> oneCellOnly = {-1, -1})

      // Member initialization list to initialize parent class (AsymmetricalJunctorInterface)
      : coupling::filtering::AsymmetricalJunctorInterface<dim>(inputCellVector1, outputCellVector1, inputCellVector2, filteredValues, "AWTF-J") {

    // Dummy write to file instance for proper junctor output. Here, only one value is written to file, i.e. _oneCellOnly = 1.
    coupling::filtering::AsymmetricalJunctorInterface<dim>::_filter1 =
        new coupling::filtering::WriteToFile<dim>(inputCellVector1, outputCellVector1, filteredValues, location[0], overwrite[0], 1);

    // Actual write to file instance for entire MD domain to file. TODO DETERMINE MD DOMAIN
    coupling::filtering::AsymmetricalJunctorInterface<dim>::_filter2 =
        new coupling::filtering::WriteToFile<dim>(inputCellVector2, {}, filteredValues, location[1], overwrite[1], oneCellOnly[1]);

  }

  ~AsymmetricalWriteToFileJunctor() {
#ifdef DEBUG_WTF_JUNCTION
    std::cout << "    AWTF-J: Destroyed AsymmetricalWriteToFileJunctor instance." << std::endl;
#endif
  }
};
