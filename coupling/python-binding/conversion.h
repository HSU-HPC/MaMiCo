// This file is part of the MaMiCo project. For conditions of distribution
// and use, please see the copyright notice in MaMiCo's main folder

#pragma once

#include <array>
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <vector>

/*
 * Functions of this namespace provide conversion bewtween popular C++ and
 * Python data structures. To be used by mamico.cpp
 * @author Felix Maurer
 */

namespace py = pybind11;

/*TODO:
 *  template for dim
 *  template for T
 */
namespace coupling {
namespace conversion {

// Converts 3-dimensional double vectors to numpy arrays
// template<class T>
py::array_t<double> stlVectorToNumpyArray_Scalar(const std::vector<double> &stl_vector, const std::vector<std::array<unsigned int, 3>> &indices) {
  if (indices.empty() || stl_vector.empty())
    throw std::runtime_error("One or more input vector is empty");

  std::vector<unsigned int> shape{(indices.back()[0]) + 1, (indices.back()[1]) + 1, (indices.back()[2]) + 1};
  py::array_t<double> res = py::array_t<double>(shape);
  auto res_unchecked = res.mutable_unchecked<3>();

  unsigned int c = 0;
  for (unsigned int i = 0; i < shape[0]; i++)
    for (unsigned int j = 0; j < shape[1]; j++)
      for (unsigned int k = 0; k < shape[2]; k++) {
        res_unchecked(i, j, k) = stl_vector[c];
        c++;
      }

  return res;
}

// template<class T>
py::array_t<double> stlVectorToNumpyArray_Vector(const std::vector<std::array<double, 3>> &stl_vector,
                                                 const std::vector<std::array<unsigned int, 3>> &indices) {
  if (indices.empty() || stl_vector.empty())
    throw std::runtime_error("One or more input vector is empty");

  std::vector<unsigned int> shape({indices.back()[0] + 1, indices.back()[1] + 1, indices.back()[2] + 1, 3});
  py::array_t<double> res = py::array_t<double>(shape);
  auto res_unchecked = res.mutable_unchecked<4>();

  unsigned int c = 0;
  for (unsigned int i = 0; i < shape[0]; i++)
    for (unsigned int j = 0; j < shape[1]; j++)
      for (unsigned int k = 0; k < shape[2]; k++) {
        res_unchecked(i, j, k, 0) = stl_vector[c][0];
        res_unchecked(i, j, k, 1) = stl_vector[c][1];
        res_unchecked(i, j, k, 2) = stl_vector[c][2];
        c++;
      }

  return res;
}

// Other way around. Still restricted to 3 dimensions and double.
// template<class T>
std::vector<double> numpyArrayToStlVector_Scalar(const py::array_t<double> &numpy_array) {
  if (numpy_array.ndim() != 3)
    throw std::runtime_error("Input array must be of exactly 3 dimensions.");

  std::vector<double> res;
  auto np_unchecked = numpy_array.unchecked<3>();

  for (unsigned int i = 0; i < numpy_array.shape(0); i++)
    for (unsigned int j = 0; j < numpy_array.shape(1); j++)
      for (unsigned int k = 0; k < numpy_array.shape(2); k++)
        res.push_back(np_unchecked(i, j, k));

  return res;
}

// template<class T>
std::vector<std::array<double, 3>> numpyArrayToStlVector_Vector(const py::array_t<double> &numpy_array) {
  if (numpy_array.ndim() != 4)
    throw std::runtime_error("Input array must be of exactly 4 dimensions.");
  if (numpy_array.shape(3) != 3)
    throw std::runtime_error("Input array's 4th dimension must be {0..2}.");

  std::vector<std::array<double, 3>> res;
  auto np_unchecked = numpy_array.unchecked<4>();

  for (unsigned int i = 0; i < numpy_array.shape(0); i++)
    for (unsigned int j = 0; j < numpy_array.shape(1); j++)
      for (unsigned int k = 0; k < numpy_array.shape(2); k++) {
        std::array<double, 3> vec = {np_unchecked(i, j, k, 0), np_unchecked(i, j, k, 1), np_unchecked(i, j, k, 2)};
        res.push_back(vec);
      }

  return res;
}

// template<class T>
std::function<std::vector<double>(std::vector<double> /*stl_vector*/, std::vector<std::array<unsigned int, 3>> /*indices*/)> *
functionWrapper_Scalar(std::function<py::array_t<double>(py::array_t<double>)> *py_func_ptr) {
  // case: py_func exists
  if (py_func_ptr) {
    // create copy at current scope
    auto py_func = *py_func_ptr;
    return new std::function<std::vector<double>(std::vector<double> /*stl_vector*/, std::vector<std::array<unsigned int, 3>> /*indices*/)>{
        [py_func](std::vector<double> stl_vector, std::vector<std::array<unsigned int, 3>> indices) {
          py::array_t<double> np_input = conversion::stlVectorToNumpyArray_Scalar(stl_vector, indices);
          py::array_t<double> np_output = py_func(np_input);
          return conversion::numpyArrayToStlVector_Scalar(np_output);
        }};
  }
  // case: py_func is "None"
  else
    return new std::function<std::vector<double>(std::vector<double> /*stl_vector*/, std::vector<std::array<unsigned int, 3>> /*indices*/)>{
        [](std::vector<double> stl_vector, std::vector<std::array<unsigned int, 3>> indices) { return stl_vector; }};
}

// template<class T>
std::function<std::vector<std::array<double, 3>>(std::vector<std::array<double, 3>> /*stl_vector*/, std::vector<std::array<unsigned int, 3>> /*indices*/)> *
functionWrapper_Vector(std::function<py::array_t<double>(py::array_t<double>)> *py_func_ptr) {
  // case: py_func exists
  if (py_func_ptr) {
    // create copy at current scope
    auto py_func = *py_func_ptr;
    return new std::function<std::vector<std::array<double, 3>>(std::vector<std::array<double, 3>> /*stl_vector*/,
                                                                std::vector<std::array<unsigned int, 3>> /*indices*/)>{
        [py_func](std::vector<std::array<double, 3>> stl_vector, std::vector<std::array<unsigned int, 3>> indices) {
          auto np_input = conversion::stlVectorToNumpyArray_Vector(stl_vector, indices);
          auto np_output = py_func(np_input);
          return conversion::numpyArrayToStlVector_Vector(np_output);
        }};
  }
  // case: py_func is "None"
  else
    return new std::function<std::vector<std::array<double, 3>>(std::vector<std::array<double, 3>> /*stl_vector*/,
                                                                std::vector<std::array<unsigned int, 3>> /*indices*/)>{
        [](std::vector<std::array<double, 3>> stl_vector, std::vector<std::array<unsigned int, 3>> indices) { return stl_vector; }};
}
} // namespace conversion
} // namespace coupling
