// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <array>

/*
 * Functions of this namespace provide conversion bewtween popular C++ and Python data structures.
 * To be used by mamico.cpp
 * @author Felix Maurer
 */

namespace coupling { namespace conversion {

	//Converts 3-dimensional double vectors to numpy arrays
	template<class T>
	py::array_t<T>* stlVectorToNumpyArray_Scalar(std::vector<T> stl_vector, std::vector<std::array<unsigned int, 3>> indices) {
		std::vector<unsigned int> shape({indices.back()[0]+1, indices.back()[1]+1, indices.back()[2]+1});
		py::array_t<T> res(shape);

		unsigned int i,j,k = 0;
		for(auto element : stl_vector){
			res[i][j][k] = element;
			(i < shape[0]-1) ? i++ : i = 0;
			(j < shape[1]-1) ? j++ : j = 0;
			(k < shape[2]-1) ? k++ : k = 0;
		}

		return res;
	}
	
	//TODO: handle arrays: 6 dimensions?
	template<class T>
	py::array_t<T>* stlVectorToNumpyArray_Vector(std::vector<std::array<T, 3>> stl_vector);
	

	//Other way around. Still restricted to 3 dimensions.
	template<class T>
	std::vector<T> numpyArrayToStlVector_Scalar(py::array_t<T> numpy_array);
	template<class T>
	std::vector<T> numpyArrayToStlVector_Vector(py::array_t<T> numpy_array);


	template<class T>
	std::vector<double> functionWrapper_Scalar(
			std::function<py::array_t<T> (py::array_t<T>)> py_func, 
			std::vector<double> stl_vector, 
			std::vector<std::array<unsigned int, 3>> indices
	){
		auto np_input = stlVectorToNumpyArray_Scalar<double>(stl_vector, indices);
		auto np_output = py_func(np_input);
		return numpyArrayToStlVector_Scalar<double>(np_output);
	}

	template<class T>
	std::vector<std::array<double, 3>> functionWrapper_Vector(
			std::function<py::array_t<T> (py::array_t<T>)> py_func, 
			std::vector<std::array<double, 3>> stl_vector, 
			std::vector<std::array<unsigned int, 3>> indices
	); /*{
		auto np_input = stlVectorToNumpyArray_Vector<double>(cells_s, indices);
		auto np_output = py_func(np_input);
		return numpyArrayToStlVector_Vector<double>(np_output);
	}*/

}}
