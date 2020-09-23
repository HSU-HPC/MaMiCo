// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling{
	template<unsigned int dim>
	class FilterFromFunction;
}


/*
 * Extension of FilterInterface.h to allow usage of custom filters using only two apply functions. 
 * Especially ment to be used for application of filters written in Python.
 *
 * Two std::function pointers are required: One for scalar and one for vector processing.
 * @author Felix Maurer
 */

template <unsigned int dim>
class coupling::FilterFromFunction : public coupling::FilterInterface<dim>{
	public:
		FilterFromFunction(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
					std::array<bool, 7> filteredValues,
					std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, dim>>)> applyScalar,
					std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)> applyVector
					):
			coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "FFF"),
			_applyScalar(applyScalar),
			_applyVector(applyVector)
		{}

		~FilterFromFunction(){}

		void operator()(){	
			
			std::vector<double> input_s;
			std::vector<std::array<double, dim>> input_v;
			std::vector<std::array<unsigned int, dim>> indices;

			//SCALAR
			for(unsigned int s = 0; s < coupling::FilterInterface<dim>::_scalarSetters.size(); s++) {
				//PACK
				for(auto cell : coupling::FilterInterface<dim>::_inputCells){
					input_s.push_back((cell->*(coupling::FilterInterface<dim>::_scalarGetters[s]))());
				}

				//APPLY
				std::vector<double> output_s = _applyScalar(input_s, indices);

				//UNPACK
				for(unsigned int i = 0; i < coupling::FilterInterface<dim>::_inputCells.size(); i++){
					(coupling::FilterInterface<dim>::_outputCells[i]->*(coupling::FilterInterface<dim>::_scalarSetters[s]))(output_s[i]);
				}	
			}

			//VECTOR
			for(unsigned int v = 0; v < coupling::FilterInterface<dim>::_vectorSetters.size(); v++) {
				//PACK
				for(auto cell : coupling::FilterInterface<dim>::_inputCells){
					tarch::la::Vector<dim, double> mamico_vec = (cell->*(coupling::FilterInterface<dim>::_vectorGetters[v]))();
					std::array<double, dim> array_vec;
					for(unsigned int d = 0; d < dim; d++) array_vec[d] = mamico_vec[d];
					input_v.push_back(array_vec);
				}

				//APPLY
				std::vector<std::array<double, dim>> output_v = _applyVector(input_v, indices);

				//UNPACK
				for(unsigned int i = 0; i < coupling::FilterInterface<dim>::_inputCells.size(); i++){
					tarch::la::Vector<dim, double> mamico_vec(-1.0);
					for(unsigned int d = 0; d < dim; d++) mamico_vec[d] = output_v[i][d];
					(coupling::FilterInterface<dim>::_outputCells[i]->*(coupling::FilterInterface<dim>::_vectorSetters[v]))(mamico_vec);
				}
			}
		}

	private:
		//this encodes what filter to use
		std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, dim>>)> _applyScalar;
		std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)> _applyVector;
};
