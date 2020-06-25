// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/IndexConversion.h"


namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}


template<unsigned int dim>
class coupling::FilterInterface{
	public:
		virtual ~FilterInterface(){/*TODO*/};
		virtual void apply(
				const std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& inputCellVector,
				std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& outputCellVector,
			   	const std::vector<tarch::la::Vector<dim,unsigned int>> cellIndices) = 0;	
};

