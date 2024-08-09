// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_ASYM_JUNCTOR_INTERFACE

#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
namespace filtering {
template <class Container_T, unsigned int dim> class AsymmetricalJunctorInterface;
}
} // namespace coupling

/**
 * Interface for Junctors in Asymmetrical Filter Junctions.
 * Cf. filtering/sequencing/AsymmetricalFilterJunction.h
 *
 * Implemenents FI. The underlying FI has the junctor's first cell data set.
 * @author Felix Maurer
 */
template <class Container_T, unsigned int dim> class coupling::filtering::AsymmetricalJunctorInterface : public coupling::filtering::FilterInterface<Container_T, dim> {
public:
  AsymmetricalJunctorInterface(
      // first cell data set
      const Container_T inputCellVector1,
      const Container_T outputCellVector1,
      const coupling::datastructures::FlexibleCellContainer<dim> inputCellVector2,
      // no output

      // parameters not specific to either 1 or 2
      const std::array<bool, 7> filteredValues, const char* type)
      : // The first cell data set in stored in FI's member variables...
        coupling::filtering::FilterInterface<Container_T, dim>(inputCellVector1, outputCellVector1, filteredValues, type),
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
  void updateCellData(Container_T& new_inputCellVector1,
                      Container_T& new_outputCellVector1,
                      coupling::datastructures::FlexibleCellContainer<dim>& new_inputCellVector2) {
    std::cout << "		AJI: Updating cell data." << std::endl;
    _inputCellVector2 = new_inputCellVector2;

    // Assumes the input c-style vectors to be nonempty. May be problematic.
    coupling::filtering::FilterInterface<Container_T, dim>::updateCellData(new_inputCellVector1, new_outputCellVector1);
  }

protected:
  /**
   * The first data set is stored in FilterInstance, the second one in here.
   * Note that the second set contains no output vector. Confer interface
   * comment above.
   */
  coupling::datastructures::FlexibleCellContainer<dim> _inputCellVector2;

  /*
   * The first cell data set should be fed to _filter1 and the second one to
   * _filter2. Note that you have to do this manually in your implementation of
   * this interface. For an example, cf. filters/WriteToFileJunctor.h
   */
  coupling::filtering::FilterInterface<Container_T, dim>* _filter1;
  coupling::filtering::FilterInterface<Container_T, dim>* _filter2;
};
