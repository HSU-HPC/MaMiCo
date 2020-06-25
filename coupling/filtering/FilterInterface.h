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

		//Applies the filter to all parts of outputCellVector() that are within the filter's domain.
		virtual void apply(
				const std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& inputCellVector, //TODO: remove this?
				std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& outputCellVector,
			   	const std::vector<tarch::la::Vector<dim,unsigned int>> cellIndices) = 0;

		//Converts indices, initializes domain_* member variables. If a filter always uses the entire MD domain, you do not need to implement this.
		//Called from coupling::FilterSequence<dim>::fillSequenceData(...)
		//virtual void initDomain();

		//Wecause we don't want to call initDomain() each time coupling::FilterSequence<dim>::fillSequenceData(...) is called,
		//we need an option to check if it's been loaded before.
		bool isDomainLoaded() { return _domainLoaded;}
	private:
		//filter domain is determined by spanning a <dim>-dimensional space between two vectors.
		//initiliazed within the filter's constructor, uses XML attributes of the corresponding filter node as input
		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;

		
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCells;
		std::vector<tarch::la::Vector<dim,unsigned int>> _domainCellIndices;

		bool _domainLoaded;
};

