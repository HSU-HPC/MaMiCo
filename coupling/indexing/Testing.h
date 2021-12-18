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
	template<unsigned int dim, coupling::indexing::IndexTrait ...>
	void checkTrivialConversions();
}

//include implementation
#include "Testing.cpph"
