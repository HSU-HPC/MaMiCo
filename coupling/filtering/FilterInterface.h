// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER

namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}

/**
 *  Generic interface for filters that are to be applied to data of coupling::MacroscopicCells before MD to Macro transfer.
 *  Examples for such filters can be found in coupling/filtering/filters.
 *  @Author Felix Maurer
 */
template<unsigned int dim>
class coupling::FilterInterface{
	public:

		/*
		 * Filter constructors are called during instanciation of their corresponding FilterSequence.
		 * You can customize parameterization in coupling::FilterSequence::loadFiltersFromXML(...).
		 */
		FilterInterface(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices):
				
				_inputCells(inputCellVector),
				_outputCells(outputCellVector),
				_cellIndices(cellIndices)
		{}


		virtual ~FilterInterface(){};

		
		//Applies the filter to all cells that are within the filter's sequence's domain.
		virtual void operator()() = 0;
	protected:
		/**
		 *  Filters should read from input vector and write to output vector.
		 *  Both vectors use the same indexing by default. 
		 *	All unmodified cells of the output vector are implicitly copied from their respective input counterpart,
		 *	i.e it is not mandatory to have any output.
		 */
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells; 	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells;
		std::vector<tarch::la::Vector<dim,unsigned int>> _cellIndices;
};
