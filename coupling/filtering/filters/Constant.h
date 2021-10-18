// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling{
	template<unsigned int dim>
	class Constant;
}

/*
 * Filter applying a constant floating point value to all cells for all filtered values.
 *
 * @author Felix Maurer
 */

template <unsigned int dim>
class coupling::Constant : public coupling::FilterInterface<dim> {
	public:
		Constant(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
					const std::array<bool, 7> filteredValues,
					const double constant
		):
			coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "Constant"),
			_constant(constant)
		{}

		void operator()() {
			for(auto outputCell : FilterInterface<dim>::_outputCells) {
				//apply to scalars
				for(auto scalarSetter : FilterInterface<dim>::_scalarSetters) {
					(outputCell->*scalarSetter)(_constant);
				}

				//apply to vectors	
				for(auto vectorSetter : FilterInterface<dim>::_vectorSetters) {
					(outputCell->*vectorSetter)({ _constant });
				}
			}
		}
	private:
		const double _constant;
};
