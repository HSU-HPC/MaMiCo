// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <cmath>
#include <string>
#include <vector>

//#define DEBUG_GAUSS
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class Gauss;

// cf. member variable in coupling::Gauss for more details
enum GaussExtrapolationStrategy { NONE, MIRROR, REFLECT };
} // namespace filtering
} // namespace coupling

// Define kernel radius. e.g. radius = 1 means kernel size of 3
#define GAUSS_KERNEL_RADIUS 1

/*
 * Gaussian filter.
 * Operates in one dimension: If you wish to use a multidimensional gaussian
 * filter, simply chain multiple instances of this filter in one FilterSequence.
 *
 * @author Felix Maurer
 */
template <unsigned int dim>
class coupling::filtering::Gauss
    : public coupling::filtering::FilterInterface<dim> {
  using coupling::filtering::FilterInterface<dim>::_inputCells;
  using coupling::filtering::FilterInterface<dim>::_outputCells;
  using coupling::filtering::FilterInterface<dim>::_scalarAccessFunctionPairs;
  using coupling::filtering::FilterInterface<dim>::_vectorAccessFunctionPairs;

  using ScalarIndex =
      coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::local,
                                    coupling::indexing::IndexTrait::md2macro,
                                    coupling::indexing::IndexTrait::noGhost>;
  using VectorIndex =
      coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector,
                                    coupling::indexing::IndexTrait::local,
                                    coupling::indexing::IndexTrait::md2macro,
                                    coupling::indexing::IndexTrait::noGhost>;

public:
  Gauss(const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
            &inputCellVector,
        const std::vector<coupling::datastructures::MacroscopicCell<dim> *>
            &outputCellVector,
        const std::array<bool, 7> filteredValues, unsigned int dimension,
        int sigma, const char *extrapolationStrategy)
      : coupling::filtering::FilterInterface<dim>(
            inputCellVector, outputCellVector, filteredValues, "GAUSS"),
        _dim(dimension), _sigma(sigma), _kernel(generateKernel()) {
    // TODO
    if (GAUSS_KERNEL_RADIUS != 1)
      throw std::runtime_error(
          "ERROR: GAUSS: Kernel radius != 1 currently not supported.");

    if (extrapolationStrategy == nullptr ||
        std::strcmp(extrapolationStrategy, "none") == 0)
      _extrapolationStrategy = NONE;
    else if (std::strcmp(extrapolationStrategy, "mirror") == 0)
      _extrapolationStrategy = MIRROR;
    else if (std::strcmp(extrapolationStrategy, "reflect") == 0)
      _extrapolationStrategy = REFLECT;
    else {
      std::cout << "Extrapolation strategy: " << extrapolationStrategy
                << std::endl;
      throw std::runtime_error("ERROR: GAUSS: Unknown extrapolation strategy.");
    }

#ifdef DEBUG_GAUSS
    std::cout << "		GAUSS (Dim: " << _dim
              << "): Created Gaussian filter." << std::endl;
    if (_extrapolationStrategy == NONE)
      std::cout << "		It will not use extrapolation." << std::endl;
    if (_extrapolationStrategy == MIRROR)
      std::cout << "		It will use mirroring extrapolation."
                << std::endl;
    if (_extrapolationStrategy == REFLECT)
      std::cout << "		It will use reflecting extrapolation."
                << std::endl;
#endif
  }

  ~Gauss() {
#ifdef DEBUG_GAUSS
    std::cout << "		GAUSS (Dim: " << _dim
              << "): Gaussian filter deconstructed" << std::endl;
#endif
  }

  void operator()();

  using CellIndex_T =
      coupling::indexing::CellIndex<dim, coupling::indexing::IndexTrait::vector,
                                    coupling::indexing::IndexTrait::local,
                                    coupling::indexing::IndexTrait::md2macro>;

private:
  std::array<double, 1 + 2 * GAUSS_KERNEL_RADIUS> generateKernel();

  constexpr double gaussianDensityFunction(int x);

  /*
   * Returns the index of the cell cell that's above the cell at index on the
   * d-axis If no such index exists, index (the first parameter) is returned.
   *
   * Index is assumed to be in terms of the MD2Macro domain, i.e. (0,..0) is the
   * lowest cell that gets sent from MD to Macro.
   */

  VectorIndex getIndexAbove(const VectorIndex index, unsigned int d);

  /*
   * Returns the index of the cell that's below the cell at index on the d-axis
   * If no such index exists, index (the first parameter) is returned.
   *
   * Index is assumed to be in terms of the MD2Macro domain, i.e. (0,..0) is the
   * lowest cell that gets sent from MD to Macro.
   */

  VectorIndex getIndexBelow(const VectorIndex index, unsigned int d);

  // on which axis this filter operates. 0 <= _dim <= dim
  const unsigned int _dim;

  // standard deviation used
  const double _sigma;

  std::array<double, 1 + 2 * GAUSS_KERNEL_RADIUS> _kernel;

  /**
   * Determines how to apply filter to border cells:
   * NONE = only use existing cells and normalize their weight accordingly
   * MIRROR = b | a b c d | c
   * REFLECT = a | a b c d | d
   *
   * The last two are congruent to SciPy's gaussian filter's respective
   * extrapolation modes
   */
  coupling::filtering::GaussExtrapolationStrategy _extrapolationStrategy;
};

// include implementation of header
#include "Gauss.cpph"
