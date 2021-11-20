// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <vector>

#include "coupling/indexing/CellIndex.h"

//#define DEBUG_FILTER_INTERFACE

namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}

using coupling::indexing::CellIndex;
using coupling::indexing::IndexTrait;

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
				const std::array<bool, 7> filteredValues,
				const char* type):
				
				_inputCells(inputCellVector),
				_outputCells(outputCellVector),
				_type(type)
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
		}

		FilterInterface(const char* type) : _type(type) {/* Used by incomplete implementations of FilterInterface. Should be redesigned via meta class.*/}


		virtual ~FilterInterface(){};

		

		/*
		 * Applies the filter to all cells that are within the filter's sequence's domain.
		 *
		 * It is very important that this method provides complete output data,
		 * i.e uses all elements of _scalarSetters and _vectorSetters on all elements of _outputCells.
		 * If this is not the case, you dont want to use this interface, but rather coupling::FilterInterfaceReadOnly
		 * and use its method copyInputToOutput().
		 */
		virtual void operator()() = 0;

		void updateCellData(
			const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& new_inputCells,
			const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& new_outputCells
		) {
			if(new_inputCells.size() != new_outputCells.size())
			   	throw std::runtime_error("New input-, output-, and indexing vectors must be of identical size.");

			_inputCells = new_inputCells;
			_outputCells = new_outputCells;
		
			#ifdef DEBUG_FILTER_INTERFACE
			std::cout << "		FI: Updated cell data." << std::endl;
			#endif
		}

		/*
		 * Basic Getters/Setters
		 */
		const char* getType() const { return _type; }
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > getInputCells() const { return _inputCells; }
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > getOutputCells() const { return _outputCells; }

		/*
		 * Advanced Getters/Setters
		 */
		coupling::datastructures::MacroscopicCell<dim>* getInputCellOfIndex(const CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>& index) {
			if(index.get() < _inputCells.size) {
				return _inputCells[index.get()];
			}
			else {
				std::cout << "Index not found: " << index << std::endl;
				throw std::runtime_error("FilterInterface: getInputCellofIndex(): Could not find index.");
			}
		}
		coupling::datastructures::MacroscopicCell<dim>* getOutputCellOfIndex(const CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>& index) {
			if(index.get() < _outputCells.size) {
				return _outputCells[index.get()];
			}
			else {
				std::cout << "Index not found: " << index << std::endl;
				throw std::runtime_error("FilterInterface: getOutputCellofIndex(): Could not find index.");
			}

		}

		/*
		 * Only used in one scenario:
		 *  - this is at index 0 in FS
		 *  - new filter gets dynamically linked into FS at index 0
		 * In that case, this was previously getting input from MD but won't any longer.
		 * The newly added filter will provide input for this one instead.
		 */
		void setInputCells(const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& newInputCells) { _inputCells = newInputCells; }

		//Size = number of cells in this filter.
		int getSize() const { return _inputCells.size(); }
		
	protected:
		/**
		 *  Filters should read from input vector and write to output vector.
		 *  Both vectors use the same indexing by default. 
		 *	All unmodified cells of the output vector are implicitly copied from their respective input counterpart,
		 *	i.e it is not mandatory to have any output.
		 */
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells;
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells;

		//scalars
		std::vector<void (coupling::datastructures::MacroscopicCell<dim>::*)(const double&)> _scalarSetters;
		std::vector<const double& (coupling::datastructures::MacroscopicCell<dim>::*)() const> _scalarGetters;
		//vectors
		std::vector<void (coupling::datastructures::MacroscopicCell<dim>::*)(const tarch::la::Vector<dim, double>&)> _vectorSetters;
		std::vector<const tarch::la::Vector<dim, double>& (coupling::datastructures::MacroscopicCell<dim>::*)() const> _vectorGetters;

		//unique identifier per filter class
		const char* _type;
};
