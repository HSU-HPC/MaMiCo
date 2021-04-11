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
				const std::array<bool, 7> filteredValues,
				const char* type):
				
				_inputCells(inputCellVector),
				_outputCells(outputCellVector),
				_cellIndices(cellIndices),
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

			//std::cout << "		First memory adresses (I/O): " << inputCellVector[0] << " "<< outputCellVector[0] << std::endl;
			
			/*std::cout << "CONSTRUCTOR: printing all cells again:" << std::endl;
			for(auto cell : _inputCells) std::cout << cell << std::endl;
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > inputCells2 = _inputCells;
			std::cout << "post copy" << std::endl;

			std::cout << "original capacity: " << _inputCells.capacity() << std::endl;
			std::cout << "copy capacity: " << inputCells2.capacity() << std::endl;

			std::cout << "original data: " << _inputCells.data() << std::endl;
			std::cout << "copy data: " << inputCells2.data() << std::endl;

			inputCells2.reserve(216);
			std::cout << "post reserve copy" << std::endl;*/
			_inputCells.reserve(216);
			/*std::cout << "post reserve original" << std::endl;*/

			//exit(0);

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

		//TODO: changed this to pass by reference. not sure if thats safe
		void updateCellData(
			std::vector<coupling::datastructures::MacroscopicCell<dim>* >& new_inputCells,
			std::vector<coupling::datastructures::MacroscopicCell<dim>* >& new_outputCells,
			std::vector<tarch::la::Vector<dim,unsigned int>>& new_cellIndices
		) {
			//std::cout << "old cell vector size: " << _inputCells.size() << std::endl;
			//std::cout << "new cell vector size: " << new_inputCells.size() << std::endl;

			if(new_inputCells.size() != new_outputCells.size() || new_outputCells.size() != new_cellIndices.size())
			   	throw std::runtime_error("New input-, output-, and indexing vectors must be of identical size.");
			
			//Alternative 1:
			//This resizing assumes the "old" _inputCells, _outputCells and _cellIndices vectors to be of identical size.
			/*if(new_inputCells.size() != _inputCells.size()) {
				std::cout << "		FI: Detected change in cell vector size. ";

				//Alternative 1.1: This would be more desireable, but segfaults. TODO
				
				_inputCells.resize(new_inputCells.size(), nullptr);
				_outputCells.resize(new_outputCells.size(), nullptr);
				_cellIndices.resize(new_cellIndices.size(), tarch::la::Vector<dim,unsigned int>(0));
				

				//Alternative 1.2: This is slightly worse imo, but also segfaults.
				int size_diff = new_inputCells.size() - _inputCells.size();
				if(size_diff > 0) { 
					for(unsigned int i = 0; i < std::abs(size_diff); i++) {
						//append dummy value size_diff times
						_inputCells.push_back(nullptr);
						_outputCells.push_back(nullptr);
						_cellIndices.push_back(tarch::la::Vector<dim,unsigned int>(0));
					}
				}
				else {
					for(unsigned int i = 0; i < std::abs(size_diff); i++) {
						_inputCells.pop_back();
						_outputCells.pop_back();
						_cellIndices.pop_back();
					}
				}

				std::cout << "Now filters " << _inputCells.size() << " cells." << std::endl;

			}

			_inputCells = new_inputCells;
			_outputCells = new_outputCells;
			_cellIndices = new_cellIndices;*/

			//Alternative 2: Unclean alternative solution: TODO: Fix above approach
			//std::cout << "all cells:" << std::endl;
			//for(auto cell : _inputCells) std::cout << cell << std::endl;

			_inputCells.clear();
			_outputCells.clear();
			_cellIndices.clear();

			/*std::vector<coupling::datastructures::MacroscopicCell<dim>* > inputCells2 = _inputCells;
			std::cout << "post copy" << std::endl;

			std::cout << "original capacity: " << _inputCells.capacity() << std::endl;
			std::cout << "copy capacity: " << inputCells2.capacity() << std::endl;

			std::cout << "original data: " << _inputCells.data() << std::endl;
			std::cout << "copy data: " << inputCells2.data() << std::endl;

			_inputCells.reserve(216);
			std::cout << "post reserve original" << std::endl;
			inputCells2.reserve(216);
			std::cout << "post reserve copy" << std::endl;*/
			
			int i = 0;
			for(auto newicell : new_inputCells) { 
				//std::cout << newicell << " I: BEFORE: " << i << std::endl; 
				_inputCells.push_back(newicell);
				//std::cout << newicell << " I: AFTER: " << i << std::endl;
			   	i++;
		   	}
			std::cout << "POST ICELLS" << std::endl;

			i = 0;
			for(auto newocell : new_outputCells) { 
				//std::cout << newocell << " O: BEFORE: " << i << std::endl; 
				_outputCells.push_back(newocell);
				//std::cout << newocell << " O: AFTER: " << i << std::endl;
			   	i++;
		   	}
			std::cout << "POST OCELLS" << std::endl;


			//for(auto newicell : new_inputCells) _inputCells.push_back(newicell);
			//for(auto newocell : new_outputCells) _outputCells.push_back(newocell);
			for(auto newindex : new_cellIndices) _cellIndices.push_back(newindex);

			std::cout << "		FI: Updated cell data." << std::endl;
		}

		void DEBUG_RESERVE(unsigned int len) {
			std::cout << "DEBUG_RESERVE (len=" << len << "): PRE RESERVE..." << std::endl;
			_inputCells.reserve(len); 
			std::cout << "POST INPUT RESERVE." << std::endl;
			_outputCells.reserve(len); 
			std::cout << "POST OUTPUT RESERVE." << std::endl;
			_cellIndices.reserve(len); 
			std::cout << "POST INDICES RESERVE." << std::endl;
		}

		const char* getType() const { return _type; }

		std::vector<coupling::datastructures::MacroscopicCell<dim>* > getInputCells() const { return _inputCells; }
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > getOutputCells() const { return _outputCells; }
		std::vector<tarch::la::Vector<dim,unsigned int>> getCellIndices() const { return _cellIndices; }

		/*
		 * Only used in one scenario:
		 *  - this is at index 0 in FS
		 *  - new filter gets dynamically linked into FS at index 0
		 * In that case, this was previously getting input from MD but won't any longer.
		 * The newly added filter will provide input for this one instead.
		 */
		void setInputCells(const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& newInputCells) { _inputCells = newInputCells; }

		//Size = number of cells in this filter.
		int getSize() const { return _cellIndices.size(); }
		
	protected:
		void DEBUG_PRINT_CELL_VELOCITY(const char* caller, unsigned int index = 0) {
			std::cout << "		" << caller << " IN ("<< (_inputCells[index]) <<"): " << _inputCells[index]->getCurrentVelocity() << std::endl;
			std::cout << "		" << caller << " OUT ("<< (_outputCells[index]) <<"): " << _outputCells[index]->getCurrentVelocity() << std::endl;
		}

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

		//unique identifier per filter class
		const char* _type;
};
