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

	template<unsigned int dim, coupling::indexing::IndexTrait ... traits, coupling::indexing::IndexTrait... convert_to_traits>
	void checkTrivialConversion(CellIndex<dim, traits...>);

}

//include implementation
#include "Testing.cpph"
