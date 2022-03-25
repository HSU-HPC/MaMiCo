// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_ASYM_JUNCTOR_INTERFACE

#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class AsymmetricalJunctorInterface;
}
} // namespace coupling

/**
 * Interface for Junctors in Asymmetrical Filter Junctions.
 * Cf. filtering/sequencing/AsymmetricalFilterJunction.h
 *
 * Implemenents FI. The underlying FI has the junctor's first cell data set.
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::filtering::AsymmetricalJunctorInterface : public coupling::filtering::FilterInterface<dim> {
public:
  AsymmetricalJunctorInterface(
      // first cell data set
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector1,
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector1,

      // first cell data set
      const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector2,
      // no output

      // parameters not specific to either 1 or 2
      const std::array<bool, 7> filteredValues, const char *type)
      : // The first cell data set in stored in FI's member variables...
        coupling::filtering::FilterInterface<dim>(inputCellVector1, outputCellVector1, filteredValues, type),
        //...while the second cell data set is stored within this class
        _inputCellVector2(inputCellVector2) {}

  virtual void operator()() {
    (*_filter1)();
    (*_filter2)();
  }

  ~AsymmetricalJunctorInterface() {
    // TODO: check for corrupt/memory leaks here
    delete _filter1;
    delete _filter2;
  };

  // TODO: make FI's updateCellData function unusable
  void updateCellData(std::vector<coupling::datastructures::MacroscopicCell<dim> *> &new_inputCellVector1,
                      std::vector<coupling::datastructures::MacroscopicCell<dim> *> &new_outputCellVector1,
                      std::vector<coupling::datastructures::MacroscopicCell<dim> *> &new_inputCellVector2) {
    std::cout << "		AJI: Updating cell data." << std::endl;
    _inputCellVector2 = new_inputCellVector2;

    // Assumes the input c-style vectors to be nonempty. May be problematic.
    coupling::filtering::FilterInterface<dim>::updateCellData(new_inputCellVector1, new_outputCellVector1);
  }

protected:
  /**
   * The first data set is stored in FilterInstance, the second one in here.
   * Note that the second set contains no output vector. Confer interface
   * comment above.
   */
  std::vector<coupling::datastructures::MacroscopicCell<dim> *> _inputCellVector2;

  /*
   * The first cell data set should be fed to _filter1 and the second one to
   * _filter2. Note that you have to do this manually in your implementation of
   * this interface. For an example, cf. filters/WriteToFileJunctor.h
   */
  coupling::filtering::FilterInterface<dim> *_filter1;
  coupling::filtering::FilterInterface<dim> *_filter2;
};
