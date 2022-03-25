// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_LA_SCALAROPERATIONS_H_
#define _TARCH_LA_SCALAROPERATIONS_H_
#include <cmath>
/**	namespace tarch */
namespace tarch {
/** namespace la */
namespace la {

/**	This function template checks if the difference of two given numbers is
 *smaller than a given tolerance i.e. if they are mathematically equal.
 *  @tparam T data type
 *  @param l First number
 *  @param r Second number
 *  @param tolerance Tolerance
 *	@return 1: if they are equal. 0: if the are not equal.
 *  @author Philipp Neumann
 */

template <class T> inline bool equals(const T &l, const T &r, const T &tolerance) { return std::abs(l - r) <= tolerance; }

} // namespace la
} // namespace tarch
#endif
