// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <array>
#include <vector>

#include "coupling/indexing/CellIndex.h"

// #define DEBUG_FILTER_INTERFACE

namespace coupling {
namespace filtering {
template <class Container_T, unsigned int dim> class FilterInterface;
}
} // namespace coupling

/**
 * Generic interface for filters that are to be applied to
 * coupling::CellContainer before MD to Macro transfer. Implementations can
 * be found in coupling/filtering/filters.
 *
 * If you wish to use a filter that does not alter its input data, i.e that is
 * read-only, you want to use coupling::FilterInterfaceReadOnly<dim> instead (as
 * provided in header file coupling/filtering/FilterPipelineReadOnly.h).
 * Examples for such filters are WriteToFile or Strouhal (in
 * coupling/filtering/filters).
 *
 * @author Felix Maurer
 */
template <class Container_T, unsigned int dim> class coupling::filtering::FilterInterface {
public:
  /*
   * Filter constructors are called during instanciation of their corresponding
   * FilterSequence. You can customize parameterization in
   * coupling::FilterSequence::loadFiltersFromXML(...).
   */
  FilterInterface(const Container_T& inputCellVector,
                  const Container_T& outputCellVector, const std::array<bool, 7> filteredValues, const char* type)
      :

        _inputCells(inputCellVector), _outputCells(outputCellVector), _type(type) {
    // microscopic mass
    if (filteredValues[0]) {
      _scalarAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getMicroscopicMass, &coupling::datastructures::CouplingCell<dim>::setMicroscopicMass});
    }
    // microscopic momentum
    if (filteredValues[1]) {
      _vectorAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getMicroscopicMomentum, &coupling::datastructures::CouplingCell<dim>::setMicroscopicMomentum});
    }
    // macroscopic mass
    if (filteredValues[2]) {
      _scalarAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getMacroscopicMass, &coupling::datastructures::CouplingCell<dim>::setMacroscopicMass});
    }
    // macroscopic momentum
    if (filteredValues[3]) {
      _vectorAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getMacroscopicMomentum, &coupling::datastructures::CouplingCell<dim>::setMacroscopicMomentum});
    }
    // potential energy
    if (filteredValues[4]) {
      _scalarAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getPotentialEnergy, &coupling::datastructures::CouplingCell<dim>::setPotentialEnergy});
    }
    // velocity
    if (filteredValues[5]) {
      _vectorAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getCurrentVelocity, &coupling::datastructures::CouplingCell<dim>::setCurrentVelocity});
    }
    // temperature
    if (filteredValues[6]) {
      _scalarAccessFunctionPairs.push_back(
          {&coupling::datastructures::CouplingCell<dim>::getTemperature, &coupling::datastructures::CouplingCell<dim>::setTemperature});
    }
  }

  FilterInterface(const char* type)
      : _type(type) { /* Used by incomplete implementations of FilterInterface.
                         Should be redesigned via meta class.*/
  }

  virtual ~FilterInterface() {};

  /*
   * Applies the filter to all cells that are within the filter's sequence's
   * domain.
   *
   * It is very important that this method provides complete output data,
   * i.e uses all elements of _scalarSetters and _vectorSetters on all elements
   * of _outputCells. If this is not the case, you dont want to use this
   * interface, but rather coupling::FilterInterfaceReadOnly and use its method
   * copyInputToOutput().
   */
  virtual void operator()() = 0;

  void updateCellData(const Container_T& new_inputCells,
                      const Container_T& new_outputCells) {
    if (new_inputCells.size() != new_outputCells.size())
      throw std::runtime_error("New input-, output-, and indexing vectors must "
                               "be of identical size.");

    _inputCells = new_inputCells;
    _outputCells = new_outputCells;

#ifdef DEBUG_FILTER_INTERFACE
    std::cout << "		FI: Updated cell data." << std::endl;
#endif
  }

  /*
   * Basic Getters/Setters
   */
  const char* getType() const { return _type; }
  Container_T getInputCells() const { return _inputCells; }
  Container_T getOutputCells() const { return _outputCells; }

  /*
   * Advanced Getters/Setters
   */
  coupling::datastructures::CouplingCell<dim>* getInputCellOfIndex(const typename Container_T::_CellIndex_T& index) {
    return _inputCells[index];
  }
  coupling::datastructures::CouplingCell<dim>* getOutputCellOfIndex(const typename Container_T::_CellIndex_T& index) {
    return _outputCells[index];
  }

  /*
   * Only used in one scenario:
   *  - this is at index 0 in FS
   *  - new filter gets dynamically linked into FS at index 0
   * In that case, this was previously getting input from MD but won't be any
   * longer. The newly added filter will provide input for this one instead.
   */
  void setInputCells(const Container_T& newInputCells) { _inputCells = newInputCells; }

  // Size = number of cells in this filter.
  int getSize() const { return _inputCells.size(); }

  /*
   * Used by filter implementations to iterate over physical properties stored
   * in a CouplingCell.
   *
   * Examplary usage:
   * 'for(auto scalar : _scalarAccessFunctionPairs)' loops over all scalar
   * properties (e.g. macro/micro mass, temperate) filtered by this specific
   * filter.
   */
  struct ScalarAccessFunctionPair {
    const double& (coupling::datastructures::CouplingCell<dim>::*get)() const; // getter function pointer
    void (coupling::datastructures::CouplingCell<dim>::*set)(const double&);   // setter function pointer
  };
  struct VectorAccessFunctionPair {
    const tarch::la::Vector<dim, double>& (coupling::datastructures::CouplingCell<dim>::*get)() const;
    void (coupling::datastructures::CouplingCell<dim>::*set)(const tarch::la::Vector<dim, double>&);
  };

protected:
  /**
   *  Filters should read from input vector and write to output vector.
   *  Both vectors use the same indexing by default.
   *	All unmodified cells of the output vector are implicitly copied from
   *their respective input counterpart, i.e it is not mandatory to have any
   *output.
   */
  Container_T _inputCells;
  Container_T _outputCells;

  // scalars getters/setters
  std::vector<ScalarAccessFunctionPair> _scalarAccessFunctionPairs;

  // vectors getters/setters
  std::vector<VectorAccessFunctionPair> _vectorAccessFunctionPairs;

  // unique identifier per filter class
  const char* _type;
};
