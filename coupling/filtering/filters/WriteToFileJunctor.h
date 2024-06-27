// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_WTF_JUNCTION

#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/interfaces/AsymmetricalJunctorInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class WriteToFileJunctor;
}
} // namespace coupling

/**
 * Combines two WriteToFile objects into one asymmetrical FilterJunctor.
 * Used to output both primary and secondary cell data in an
 * AsymmetricalFilterJunction.
 *
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::filtering::WriteToFileJunctor : public coupling::filtering::AsymmetricalJunctorInterface<dim> {
public:
  WriteToFileJunctor(
      // first cell data set
      const std::vector<coupling::datastructures::CouplingCell<dim>*> inputCellVector1,
      const std::vector<coupling::datastructures::CouplingCell<dim>*> outputCellVector1,

      // second cell data set
      const std::vector<coupling::datastructures::CouplingCell<dim>*> inputCellVector2,
      // no secondary output cells

      //"global" parameters for both WriteToFile instances
      const std::array<bool, 7> filteredValues,

      // WriteToFile-specific parameters. [0] is for the first WriteToFile
      // instance and [1] for the second one respectively.
      std::array<std::string, 2> location, std::array<bool, 2> overwrite = {false}, std::array<int, 2> oneCellOnly = {-1})
      : coupling::filtering::AsymmetricalJunctorInterface<dim>(inputCellVector1, outputCellVector1, inputCellVector2, filteredValues, "WTF-J") {
    // write to file instance covering first cell data set
    coupling::filtering::AsymmetricalJunctorInterface<dim>::_filter1 =
        new coupling::filtering::WriteToFile<dim>(inputCellVector1, outputCellVector1, filteredValues, location[0], overwrite[0], oneCellOnly[0]);

    // write to file instance covering second cell data set
    coupling::filtering::AsymmetricalJunctorInterface<dim>::_filter2 =
        new coupling::filtering::WriteToFile<dim>(inputCellVector2, {}, // no output
                                                  filteredValues, location[1], overwrite[1], oneCellOnly[1]);
  }

  ~WriteToFileJunctor() {
#ifdef DEBUG_WTF_JUNCTION
    std::cout << "    WTF-J: Destroyed WriteToFileJunctor instance." << std::endl;
#endif
  }
};
