// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterfaceReadOnly.h"
#include <algorithm>
#include <vector>

#define DEBUG_STROUHAL

namespace coupling {
namespace filtering {
template <unsigned int dim> class Strouhal;
}
} // namespace coupling

/**
 * Implements a read-only filter evaluating given input velocity data to
 * approximate a Strouhal number in y-direction. As an alternative to this,
 * Numpy-compatible implementations of FFT can be linked into a FS using
 * FilterFromFunction. Fourier transformations are expected to yield more
 * accurate results.
 * @author Felix Maurer
 */
template <unsigned int dim> class coupling::filtering::Strouhal : public coupling::filtering::FilterInterfaceReadOnly<dim> {
public:
  Strouhal(const std::vector<coupling::datastructures::CouplingCell<dim>*>& inputCellVector,
           const std::vector<coupling::datastructures::CouplingCell<dim>*>& outputCellVector,
           const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, std::array<bool, 7> filteredValues, double u, double d)
      : coupling::filtering::FilterInterfaceReadOnly<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "STROUHALCPP"), _U(u), _D(d) {
    if (dim < 2) {
      std::cout << "ERROR: Strouhal filter only works for dim >= 2." << std::endl;
      exit(EXIT_FAILURE);
    }
#ifdef DEBUG_STROUHAL
    std::cout << "		STROUHAL: Instance created." << std::endl;
#endif
  }

  ~Strouhal() {
    std::cout << "STROUHAL NUMBER IN MD DOMAIN: " << calculateStrouhalNumber() << std::endl;
#ifdef DEBUG_WRITE_TO_FILE
    std::cout << "		STROUHAL: Instance deconstructed." << std::endl;
#endif
  }

  void operator()();

private:
  double calculateStrouhalNumber();

  std::vector<double> _v_y;
  double _U;
  double _D;
};

// include implementation of header
#include "Strouhal.cpph"

/*
 * TODO:
 * allow offset (e.g. start measuring after 500 coupling cycles)
 */
