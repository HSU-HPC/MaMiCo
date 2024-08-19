// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
#include <functional>

namespace coupling {
namespace filtering {
template <unsigned int dim> class FilterFromFunction;
}
} // namespace coupling

/*
 * Extension of FilterInterface.h to allow usage of custom filters using only
 * two apply functions. Especially ment to be used for application of filters
 * written in Python.
 *
 * Two std::function pointers are required: One for scalar and one for vector
 * processing.
 * @author Felix Maurer
 */

template <unsigned int dim> class coupling::filtering::FilterFromFunction : public coupling::filtering::FilterInterface<dim> {
public:
  FilterFromFunction(
      const std::vector<coupling::datastructures::CouplingCell<dim>*>& inputCellVector,
      const std::vector<coupling::datastructures::CouplingCell<dim>*>& outputCellVector, std::array<bool, 7> filteredValues,
      const std::function<std::vector<double>(std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* applyScalar,
      const std::function<std::vector<std::array<double, dim>>(std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* applyVector)
      : coupling::filtering::FilterInterface<dim>(inputCellVector, outputCellVector, filteredValues, "FFF"), _applyScalar(applyScalar),
        _applyVector(applyVector) {

    if (applyScalar == nullptr or applyVector == nullptr)
      throw std::runtime_error("ERROR: FilterFromFunction received nullptr as function pointer!");

    // cast MaMiCo indexing to std::array
    tarch::la::Vector<dim, int> mamicoIndex;
    std::array<unsigned int, dim> stlIndex;
    for (unsigned int i = 0; i < inputCellVector.size(); i++) {
      // interpret position of cell in inputCellVector as linear local
      // md-to-macro index, then convert it to vector
      using coupling::indexing::IndexTrait;
      mamicoIndex = coupling::indexing::convertToVector<dim>(I14{i});

      for (unsigned int d = 0; d < dim; d++)
        stlIndex[d] = mamicoIndex[d];
      _stlIndices.push_back(stlIndex);
    }
  }

  ~FilterFromFunction() {
    delete _applyScalar;
    delete _applyVector;
  }

  void operator()() {
    std::vector<double> input_s;
    std::vector<std::array<double, dim>> input_v;

    input_s.reserve(coupling::filtering::FilterInterface<dim>::_inputCells.size());
    input_v.reserve(coupling::filtering::FilterInterface<dim>::_inputCells.size());

    /*
     * SCALAR
     */
    for (const auto scalarProperty : coupling::filtering::FilterInterface<dim>::_scalarAccessFunctionPairs) {

      /*
       * PACK
       */
      for (auto cell : coupling::filtering::FilterInterface<dim>::_inputCells) {
        input_s.push_back((cell->*scalarProperty.get)());
      }

      // std::cout << "Now applying scalar func at: " << _applyScalar <<
      // std::endl;
      /*
       * APPLY
       *
       * If the cell vector is emtpy, we skip this step. It would not have any
       * effect anyway and could cause problems if the apply functions does
       * cannot handle empty input sets.
       */
      std::vector<double> output_s = {};
      if (input_s.size() > 0)
        output_s = (*_applyScalar)(input_s, _stlIndices);

      input_s.clear();

      /*
       * UNPACK
       */
      for (unsigned int i = 0; i < coupling::filtering::FilterInterface<dim>::_inputCells.size(); i++) {
        (coupling::filtering::FilterInterface<dim>::_outputCells[i]->*scalarProperty.set)(output_s[i]);
      }
    }

    /*
     * VECTOR
     */
    for (const auto vectorProperty : coupling::filtering::FilterInterface<dim>::_vectorAccessFunctionPairs) {

      // coupling::FilterInterface<dim>::DEBUG_PRINT_CELL_VELOCITY("FFF BEFORE
      // ");

      /*
       * PACK
       */
      for (auto cell : coupling::filtering::FilterInterface<dim>::_inputCells) {
        tarch::la::Vector<dim, double> mamico_vec = (cell->*vectorProperty.get)();
        std::array<double, dim> array_vec;
        for (unsigned int d = 0; d < dim; d++)
          array_vec[d] = mamico_vec[d];
        input_v.push_back(array_vec);
      }

      // std::cout << "Now applying vector func at: " << _applyVector <<
      // std::endl;
      /*
       * APPLY
       *
       * Cf. scalar case.
       */
      std::vector<std::array<double, dim>> output_v = {};
      if (input_v.size() > 0)
        output_v = (*_applyVector)(input_v, _stlIndices);

      input_v.clear();

      /*
       * UNPACK
       */
      for (unsigned int i = 0; i < coupling::filtering::FilterInterface<dim>::_inputCells.size(); i++) {
        tarch::la::Vector<dim, double> mamico_vec{};
        for (unsigned int d = 0; d < dim; d++)
          mamico_vec[d] = output_v[i][d];
        (coupling::filtering::FilterInterface<dim>::_outputCells[i]->*vectorProperty.set)(mamico_vec);
      }

      // coupling::FilterInterface<dim>::DEBUG_PRINT_CELL_VELOCITY("FFF AFTER
      // ");
    }
  }

private:
  // FFFs use slightly different datastructures for index/cell storage than
  // other filters
  std::vector<std::array<unsigned int, dim>> _stlIndices;

  // this encodes what filter to use
  const std::function<std::vector<double>(std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* _applyScalar;
  const std::function<std::vector<std::array<double, dim>>(std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* _applyVector;
};
