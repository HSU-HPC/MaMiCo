// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include <vector>

//#define DEBUG_FILTER

namespace coupling{
    template<unsigned int dim, std::size_t inputc, std::size_t outputc>
    class JunctorInterface;
}

/**
 * Implemenents FI. The underlying FI has the junctor's main partition's input/output data. (cf. 30, 55f.)
 * TODO: Further comments.
 *  @Author Felix Maurer
 */
template<unsigned int dim, std::size_t inputc, std::size_t outputc>
class coupling::JunctorInterface : public coupling::FilterInterface<dim> {
	public:
		JunctorInterface(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVectors[inputc],
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVectors[outputc],
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
				const std::array<bool, 7> filteredValues,
				const char* type):
			//This assumes the vector of cell vectors to be nonempty. Suboptimal
			coupling::FilterInterface<dim>(inputCellVectors[0], outputCellVectors[0], cellIndices, filteredValues, type),
			_inputCellVectors(inputCellVectors),
			_outputCellVectors(outputCellVectors)
		{
		}

		//TODO: do i need this constructor here as well?
		//FilterInterface(const char* type) : _type(type) {/* Used by incomplete implementations of FilterInterface. Should be redesigned via meta class.*/}


		virtual ~JunctorInterface(){};


		//TODO: changed this to pass by reference. not sure if thats safe
		void updateCellData(
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > new_inputCellVectors[inputc],
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > new_outputCellVectors[outputc],
			std::vector<tarch::la::Vector<dim,unsigned int>>& new_cellIndices
		) {
			std::cout << "		JI: Updating cell data." << std::endl;
			_inputCellVectors = new_inputCellVectors;
			_outputCellVectors = new_outputCellVectors;

			//Assumes the input c-style vectors to be nonempty. May be problematic.
			coupling::FilterInterface<dim>::updateCellData(
					new_inputCellVectors[0],
					new_outputCellVectors[0],
					new_cellIndices);
		}

	protected:
		/**
		 * Unlike regular filters, junctors allow for multiple input- and output-sets
		 */
		//TODO: is it a good idea for this to be a c-style array? cf. updateCellData. perhaps std::array instead?
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVectors[inputc]; 	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCellVectors[outputc];
};