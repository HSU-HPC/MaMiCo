// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}


template<unsigned int dim>
class coupling::FilterInterface{
	public:
		virtual std::vector<coupling::datastructures::MacroscopicCell<dim>* > apply(
				std::vector<coupling::datastructures::MacroscopicCell<dim>*  > inputCellVector,
			   	std::vector<unsigned int> cellIndices) = 0;
};

