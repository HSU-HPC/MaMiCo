// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_LA_SCALAROPERATIONS_H_
#define _TARCH_LA_SCALAROPERATIONS_H_
#include <cmath>

namespace tarch{
namespace la{

template<class T>
inline bool equals(const T& l, const T& r, const T& tolerance){
  return std::abs(l-r) <= tolerance;
}

}
}
#endif 
