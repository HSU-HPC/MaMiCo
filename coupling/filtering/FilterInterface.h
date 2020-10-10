// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <vector>

//#define DEBUG_FILTER

namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}

/**
 *  Generic interface for filters that are to be applied to data of coupling::MacroscopicCells before MD to Macro transfer.
 *  Examples for such filters can be found in coupling/filtering/filters.
 *
 *  If you wish to use a filter that does not give cell output data, i.e that is read-only, you want to use 
 *  	coupling::FilterInterfaceReadOnly<dim>
 *  instead (as provided in header file coupling/filtering/FilterPipelineReadOnly.h).
 *  Examples for such filters are WriteToFile or Strouhal (in coupling/filtering/filters).
 *
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
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
				bool filteredValues[7]):
				
				_inputCells(inputCellVector),
				_outputCells(outputCellVector),
				_cellIndices(cellIndices)
		{
			if(filteredValues[0]){
				_scalarSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setMicroscopicMass);
				_scalarGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getMicroscopicMass);
			}
			if(filteredValues[1]){
				_vectorSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setMicroscopicMomentum);
				_vectorGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getMicroscopicMomentum);
			}
			if(filteredValues[2]){
				_scalarSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setMacroscopicMass);
				_scalarGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getMacroscopicMass);
			}
			if(filteredValues[3]){
				_vectorSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setMacroscopicMomentum);
				_vectorGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getMacroscopicMomentum);
			}
			if(filteredValues[4]){
				_scalarSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setPotentialEnergy);
				_scalarGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getPotentialEnergy);
			}
			if(filteredValues[5]){
				_vectorSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setCurrentVelocity);
				_vectorGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getCurrentVelocity);
			}
			if(filteredValues[6]){
				_scalarSetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::setTemperature);
				_scalarGetters.push_back(&coupling::datastructures::MacroscopicCell<dim>::getTemperature);
			}

			//std::cout << "		First memory adresses (I/O): " << inputCellVector[0] << " "<< outputCellVector[0] << std::endl;
		}


		virtual ~FilterInterface(){};

		
		//Applies the filter to all cells that are within the filter's sequence's domain.

		//It is very important that this method provides complete output data,
		//i.e uses all elements of _scalarSetters and _vectorSetters on all elements of _outputCells.
		//If this is not the case, you dont want to use this interface, but rather
		//	coupling::FilterInterfaceReadOnly
		//and use its method copyInputToOutput().
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

		//scalars
		std::vector<void (coupling::datastructures::MacroscopicCell<dim>::*)(const double&)> _scalarSetters;
		std::vector<const double& (coupling::datastructures::MacroscopicCell<dim>::*)() const> _scalarGetters;
		//vectors
		std::vector<void (coupling::datastructures::MacroscopicCell<dim>::*)(const tarch::la::Vector<dim, double>&)> _vectorSetters;
		std::vector<const tarch::la::Vector<dim, double>& (coupling::datastructures::MacroscopicCell<dim>::*)() const> _vectorGetters;


};
