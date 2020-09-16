// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "FilterInterface.h"

namespace coupling{
	template<unsigned int dim>
	class FilterInterfaceReadOnly;
}

/*
 * Extension of FilterInterface.h for cases in which the filter itself does not produce any output data.
 * For such filters, you want to make use of copyInputToOutput() (see below) in every filter step.
 */

template <unsigned int dim>
class coupling::FilterInterfaceReadOnly : public coupling::FilterInterface<dim>{
	public:
		FilterInterfaceReadOnly(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
					const std::array<bool, 7> filteredValues):
			coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues)
		{}
	protected:
		/*
		 * Copies all filtered data from input to output. You always want to call this as part of any implementation of
		 * 		coupling::FilterInterface<dim>::operator()
		 * 	when implementing this interface, that is implementing a read-only filter (e.g WriteToFile, Storuhal)
		 * 	If you would not do that, the successors of the implementing filter in a sequence would get faulty input data.
		 */
		void copyInputToOutput(){
			for(unsigned int ci = 0; ci < coupling::FilterInterface<dim>::_outputCells.size(); ci++){
				for(unsigned int si = 0; si < coupling::FilterInterface<dim>::_scalarSetters.size(); si++){
					(coupling::FilterInterface<dim>::_outputCells[ci]->*(coupling::FilterInterface<dim>::_scalarSetters[si]))(
					(coupling::FilterInterface<dim>::_inputCells[ci]->*(coupling::FilterInterface<dim>::_scalarGetters[si]))());
				}
				for(unsigned int vi = 0; vi < coupling::FilterInterface<dim>::_vectorSetters.size(); vi++){
					(coupling::FilterInterface<dim>::_outputCells[ci]->*(coupling::FilterInterface<dim>::_vectorSetters[vi]))(
					(coupling::FilterInterface<dim>::_inputCells[ci]->*(coupling::FilterInterface<dim>::_vectorGetters[vi]))());
				}
			}

		}

};
