// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_LA_MATRIX_H_
#define _TARCH_LA_MATRIX_H_
#include "tarch/TarchDefinitions.h"

namespace tarch {
namespace la {
template <int rows, int cols, class T> class Matrix;
}
} // namespace tarch

/**	used to define 2-D matrices of arbitrary types.
 *  @tparam rows defines the number f rows
 *  @tparam cols defines the number f columns
 *  @tparam T data type
 *  @author Philipp Neumann
 */

template <int rows, int cols, class T> class tarch::la::Matrix {
private:
  T _entries[rows][cols];

public:
  /** @brief Two constructors, a default one and a constructor, which
initializes all elemnet of the matrix with entry t
    * 	@param t the initial value for all elemnts of the matrix.
        */
  Matrix() {}
  Matrix(const T &t) {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        _entries[i][j] = t;
      }
    }
  }

  /** @brief operator(i, j) gives back the the element (i, j)
   * 	@param i row index
   * 	@param j column index
   */
  T &operator()(int i, int j) {
#if (TARCH_DEBUG == TARCH_YES)
    if (i < 0 || i > rows - 1) {
      std::cout << "ERROR Matrix T& operator(): i out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (j < 0 || j > cols - 1) {
      std::cout << "ERROR Matrix T& operator(): j out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _entries[i][j];
  }
  const T &operator()(int i, int j) const {
#if (TARCH_DEBUG == TARCH_YES)
    if (i < 0 || i > rows - 1) {
      std::cout << "ERROR Matrix const T& operator(): i out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if (j < 0 || j > cols - 1) {
      std::cout << "ERROR Matrix const T& operator(): j out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _entries[i][j];
  }
};

#endif
