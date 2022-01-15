#pragma once

namespace coupling::indexing::testing {	

	/**
	 * TODO: comment
	 */
	template<unsigned int dim>
	void printAllBoundaries(std::ostream& = std::cout);

	/**
	 * TODO: comment
	 */
	template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
	void checkAllTrivialConversions();

	template<unsigned int dim, class Index, coupling::indexing::IndexTrait... convert_to_traits>
	void checkTrivialConversion(Index);

}

//include implementation
#include "Testing.cpph"
