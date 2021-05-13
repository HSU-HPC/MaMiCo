#pragma once



#include "tensorflow/core/framework/tensor.h"
#include "ConsoleOutput.h"
#include "ReadCSV.h"
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

float tensorMean(tensorflow::Tensor t);

std::vector<float> operator-=(std::vector<float>& vec1, std::vector<float> vec2);

float maxInColumn(std::vector<std::vector<float>> vec, int column);

const tensorflow::Tensor VecToTensor(std::vector<float> vec);

std::vector<std::vector<float>> TensorToVec(tensorflow::Tensor t);

//returns single 1d tensor out of 2d tensor
tensorflow::Tensor getTensorByIndex(tensorflow::Tensor t, int index);

//#include "TensorCalc.cpph"


