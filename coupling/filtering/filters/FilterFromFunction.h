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
					bool filteredValues[7],
					std::function<std::vector<double> (std::vector<double>, std::vector<std::array<double, dim>>)> applyScalar,
					std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<double, dim>>)> applyVector
					):
			coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues),
			_applyScalar(applyScalar),
			_applyVector(applyVector)
		{}

		~FilterFromFunction(){}

		//TODO: use filteredValues[] instead of exclusively macro mass/momentum 
		void operator()(){	
			//PACK
			std::vector<double> input_s;
			std::vector<std::array<double, dim>> input_v;
			std::vector<std::array<double, dim>> indices;

			for(auto cell : coupling::FilterInterface<dim>::_inputCells){
				//vector
				tarch::la::Vector<dim, double> mamico_vec = cell->getMacroscopicMomentum();
				std::array<double, dim> array_vec;
				for(unsigned int d = 0; d < dim; d++) array_vec[d] = mamico_vec[d];
				input_v.push_back(array_vec);

				//scalar
				input_s.push_back(cell->getMacroscopicMass());
			}

			//APPLY
			std::vector<double> output_s = _applyScalar(input_s, indices);
			std::vector<std::array<double, dim>> output_v = _applyVector(input_v, indices);

			if(output_s.size() != coupling::FilterInterface<dim>::_cellIndices.size() || output_v.size() != coupling::FilterInterface<dim>::_cellIndices.size())
				std::cout << "Something went very wrong at data structure conversion. :(";


			//UNPACK
			for(unsigned int i = 0; i < coupling::FilterInterface<dim>::_inputCells.size(); i++){
				//vector
				tarch::la::Vector<dim, double> mamico_vec(-1.0);
				for(unsigned int d = 0; d < dim; d++) mamico_vec[d] = output_v[i][d];
				coupling::FilterInterface<dim>::_outputCells[i]->setMacroscopicMomentum(mamico_vec);

				//scalar
				coupling::FilterInterface<dim>::_outputCells[i]->setMacroscopicMass(output_s[i]);
			}	
		}

	private:
		//this encodes what filter to use
		std::function<std::vector<double> (std::vector<double>, std::vector<std::array<double, dim>>)> _applyScalar;
		std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<double, dim>>)> _applyVector;
};
