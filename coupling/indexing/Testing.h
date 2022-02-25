#pragma once

namespace coupling::indexing::testing {

/**
 * TODO: comment
 */
template <unsigned int dim> void printAllBoundaries(std::ostream & = std::cout);

/**
 * TODO: comment
 */
template <unsigned int dim> void printAllDomains(std::ostream & = std::cout);

template <unsigned int dim, coupling::indexing::IndexTrait... traits>
void printIndexDomain(std::ostream & = std::cout);

/**
 * TODO: comment
 */
template <unsigned int dim, coupling::indexing::IndexTrait... traits>
void checkAllTrivialConversions();

template <unsigned int dim, class Index,
          coupling::indexing::IndexTrait... convert_to_traits>
void checkTrivialConversion(Index);

} // namespace coupling::indexing::testing

// include implementation
#include "Testing.cpph"
