/*
contains some functions to manipulate/convert tensorflow::Tensor and std::vector
*/
#pragma once

#include "tensorflow/core/framework/tensor.h"

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

//returns mean of all tensor entries
float tensorMean(tensorflow::Tensor t);

//subtracts 1d vectors
std::vector<float> operator-=(std::vector<float>& vec1, std::vector<float> vec2);

//returns largest float within the given column
float maxInColumn(std::vector<std::vector<float>> vec, int column);

//converts 2d vector to 2d tensor with the same entries
tensorflow::Tensor VecToTensor(std::vector<std::vector<float>> vec);

//converts 1d vector to 2d tensor with the same entries
//the tensor is 2d because predictions/training required 2d tensor
tensorflow::Tensor VecToTensor(std::vector<float> vec);

//converts 2d tensor to 2d vector
std::vector<std::vector<float>> TensorToVec(tensorflow::Tensor t);

//returns single entry out of 2d tensor
//can be used to get a single timestep out of a batch
//the returned entry is also a 2d tensor to allow it to be used for predictions
tensorflow::Tensor getTensorByIndex(tensorflow::Tensor t, int index);


