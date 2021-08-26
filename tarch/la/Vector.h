// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_LA_VECTOR_H_
#define _TARCH_LA_VECTOR_H_
#include<iostream>
#include <cmath>
#include <type_traits>
#include "tarch/TarchDefinitions.h"

namespace tarch { namespace la {
  template<int size, class T> class Vector;
} }


/** light-weight implementation of a vector class.
 *  @author Philipp Neumann
 */
template<int size, class T>
class tarch::la::Vector {
  private:
    T _entries[size];
  public:
    Vector(){}
    /** init. vector with a scalar value. */
    Vector(const T& t){
      for (int i = 0; i < size; i++){ _entries[i] = t; }
    }
    /** special constructors for 2D and 3D; left empty for general purpose vector*/
    Vector(const T& t0, const T& t1){
      static_assert(size==2, "ERROR Vector(t0,t1) only valid for 2D vectors!");
      _entries[0] = t0;
      _entries[1] = t1;
    }
    Vector(const T& t0, const T& t1, const T& t2){
      static_assert(size==3, "ERROR Vector(t0,t1,t2) only valid for 3D vectors!");
      _entries[0] = t0;
      _entries[1] = t1;
      _entries[2] = t2;
    }
    /** init. vector from vector */
    Vector(const Vector<size,T>& v){
      for (int i = 0; i < size; i++){ _entries[i] = v[i]; }
    }
    /** assigns a value to all vector entries. */
    void assign(const T& t){
      for (int i = 0; i < size; i++){ _entries[i] = t; }
    }
    /** assign vector */
    Vector<size,T>& operator=(const Vector<size,T> &v){
      for (int i = 0; i < size; i++) { _entries[i] = v[i]; }
      return *this;
    }

    /** access to vector entries; both () and [] are allowed */
    T& operator[]( int i) {
      #if (TARCH_DEBUG==TARCH_YES)
      if(i<0 || i>size-1){std::cout << "ERROR Vector T& operator[]: i out of range!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return _entries[i];
    }
    const T& operator[](int i) const {
      #if (TARCH_DEBUG==TARCH_YES)
      if (i<0 || i>size-1){std::cout << "ERROR Vector const T& operator[]: i out of range!" << std::endl; exit(EXIT_FAILURE);}
      #endif
      return _entries[i];
    }

    /** add a vector to this existing one */
    Vector<size,T>& operator+=(const Vector<size,T> & v){
      for (int i = 0; i < size; i++){ _entries[i] += v[i];}
      return *this;
    }
    /** subtracts a vector */
    Vector<size,T>& operator-=(const Vector<size,T> &v){
      for (int i = 0; i < size; i++){ _entries[i] -= v[i];}
      return *this;
    }
};


#include "VectorOperations.cpph"
#endif // _TARCH_LA_VECTOR_H_
