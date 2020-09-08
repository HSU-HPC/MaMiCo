// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "FilterInterface.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace coupling{
	template<unsigned int dim>
	class FilterInterfacePython;
}

namespace py = pybind11;

/*
 * Extension of FilterInterface.h to allow usage of filters written in Python.
 * @author Felix Maurer
 */

template <unsigned int dim>
class coupling::FilterInterfacePython : public coupling::FilterInterface<dim>{
	public:
		FilterInterfacePython(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
					bool filteredValues[7],
					const char* source,
					const char* pArgs):
			coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues),
			_source(source),
			_pArgs(pArgs)
		{}

		~FilterInterfacePython(){ /*TODO*/}

		//TODO: use filteredValues[] instead of exclusively macro mass/momentum 
		void operator()(){
	   		py::initialize_interpreter();
			py::module source = py::module::import(_source);

			py::function apply_scalar = py::reinterpret_borrow<py::function>(source.attr("apply_scalar"));
			py::function apply_vector = py::reinterpret_borrow<py::function>(source.attr("apply_vector"));
				
			//PACK
			std::vector<double> input_s;
			std::vector<std::vector<double>> input_v;
			for(auto cell : coupling::FilterInterface<dim>::_inputCells){
				//vector
				tarch::la::Vector<dim, double> mamico_vec = cell->getMacroscopicMomentum();
				std::vector<double> vec;
				for(unsigned int d = 0; d < dim; d++) vec.push_back(mamico_vec[d]); //cast mamico-vec -> std::vec?
				input_v.push_back(vec);

				//scalar
				input_s.push_back(cell->getMacroscopicMass());
			}

			//APPLY
			std::vector<double> output_s = apply_scalar(input_s, _pArgs).cast<std::vector<double>>();
			std::vector<std::vector<double>> output_v = apply_vector(input_v, _pArgs).cast<std::vector<std::vector<double>>>();
			if(output_s.size() != coupling::FilterInterface<dim>::_cellIndices.size() || output_v.size() != coupling::FilterInterface<dim>::_cellIndices.size())
				std::cout << "Something went very wrong at data structure conversion. :(";

			py::finalize_interpreter();

			//UNPACK
			for(unsigned int i = 0; i < coupling::FilterInterface<dim>::_inputCells.size(); i++){
				//vector
				tarch::la::Vector<dim, double> mamico_vec(-1.0);
				for(unsigned int d = 0; d < dim; d++) mamico_vec[d] = output_v[i][d];
				coupling::FilterInterface<dim>::_outputCells[i]->setMacroscopicMomentum(mamico_vec); //cast std::vec -> mamico-vec?

				//scalar
				coupling::FilterInterface<dim>::_outputCells[i]->setMacroscopicMass(output_s[i]);
			}	
		}

	protected:
		//this encodes what filter to use
		const char* _source;
		const char* _pArgs;

	private:
		//TODO: remove?
		void unpackNumpyArray(py::array_t<double> outputData);


};
