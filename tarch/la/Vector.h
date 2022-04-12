// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_LA_VECTOR_H_
#define _TARCH_LA_VECTOR_H_
#include "tarch/TarchDefinitions.h"
#include <cmath>
#include <iostream>
#include <type_traits>

namespace tarch {
namespace la {
template <int size, class T> class Vector;
}
} // namespace tarch

/**	light-weight implementation of a vector class.
 *  @tparam size Size
 *  @tparam T Data type
 *  @author Philipp Neumann
 */

template <int size, class T> class tarch::la::Vector {
private:
  T _entries[size];

public:
  Vector() {}
  /** @brief init. vector with a scalar value.
   * 	@param t scalar value for the initialization
   */

  Vector(const T& t) {
    for (int i = 0; i < size; i++) {
      _entries[i] = t;
    }
  }
  /** @brief special constructor for 2D left empty for general purpose vector
   * 	@param t0 scalar value for the initialization
   * 	@param t1 scalar value for the initialization
   */
  Vector(const T& t0, const T& t1) {
    static_assert(size == 2, "ERROR Vector(t0,t1) only valid for 2D vectors!");
    _entries[0] = t0;
    _entries[1] = t1;
  }
  /** @brief special constructor for 3D; left empty for general purpose vector
   * 	@param t0 scalar value for the initialization
   * 	@param t1 scalar value for the initialization
   * 	@param t2 scalar value for the initialization
   */
  Vector(const T& t0, const T& t1, const T& t2) {
    static_assert(size == 3, "ERROR Vector(t0,t1,t2) only valid for 3D vectors!");
    _entries[0] = t0;
    _entries[1] = t1;
    _entries[2] = t2;
  }
  /** @brief constructor init. vector from vector
   * 	@param v Vector for the initialization
   */
  Vector(const Vector<size, T>& v) {
    for (int i = 0; i < size; i++) {
      _entries[i] = v[i];
    }
  }
  /** @brief assigns a value to all vector entries.
   * 	@param t scalar value for the assignment
   */
  void assign(const T& t) {
    for (int i = 0; i < size; i++) {
      _entries[i] = t;
    }
  }
  /** @brief operator overloading for vector assignment
   * 	@param v Vector for the assignment
   */
  Vector<size, T>& operator=(const Vector<size, T>& v) {
    for (int i = 0; i < size; i++) {
      _entries[i] = v[i];
    }
    return *this;
  }
  /** @brief operator overloading; access to vector entries; both () and [] are
   * allowed
   * 	@param i index
   */
  T& operator[](int i) {
#if (TARCH_DEBUG == TARCH_YES)
    if (i < 0 || i > size - 1) {
      std::cout << "ERROR Vector T& operator[]: i out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _entries[i];
  }
  /** @brief operator overloading; access to vector entries; both () and [] are
   * allowed; !!! Attention: const
   * 	@param i index
   */
  const T& operator[](int i) const {
#if (TARCH_DEBUG == TARCH_YES)
    if (i < 0 || i > size - 1) {
      std::cout << "ERROR Vector const T& operator[]: i out of range!" << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
    return _entries[i];
  }
  /** @brief operator overloading; add a vector to this existing one (this)
   * 	@param v vector that has to be added
   */
  Vector<size, T>& operator+=(const Vector<size, T>& v) {
    for (int i = 0; i < size; i++) {
      _entries[i] += v[i];
    }
    return *this;
  }
  /** @brief operator overloading; subtracts a vector from the existing one
   * (this)
   * 	@param v vector that has to be subtracted
   */
  Vector<size, T>& operator-=(const Vector<size, T>& v) {
    for (int i = 0; i < size; i++) {
      _entries[i] -= v[i];
    }
    return *this;
  }

  /** cast Vector<size, T> to Vector<size, convert_to_T> */
  template <class convert_to_T> explicit operator Vector<size, convert_to_T>() const {
    Vector<size, convert_to_T> ans;
    for (unsigned int i = 0; i < size; i++) {
      ans[i] = static_cast<convert_to_T>(_entries[i]);
    }
    return ans;
  }
};

#include "VectorOperations.cpph"
#endif // _TARCH_LA_VECTOR_H_
