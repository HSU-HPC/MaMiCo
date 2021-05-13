// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_NeuralNet

#include "coupling/filtering/interfaces/JunctorInterface.h"


#include "tensorflow/cc/ops/array_ops.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/cc/framework/ops.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/cc/ops/state_ops.h"
#include "tensorflow/cc/ops/math_ops.h"
#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/TensorCalc.h"

namespace coupling {
  template<unsigned int dim>
  class NeuralNet;
}

/*auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize) {
  int32_t* tensor_ptr = tensor.flat<int32_t>().data();
  std::vector<int32_t> v(tensor_ptr, tensor_ptr + tensorSize);
  return v;
}
*/

template<unsigned int dim>
class coupling::NeuralNet : public coupling::AsymmetricalJunctorInterface<dim>{
	public:
		NeuralNet(
			//first cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector1,
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector1,
			const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices1,
			const std::vector<tarch::la::Vector<dim, unsigned int>> sequenceCellIndices1,

			//second cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector2,
			//no output cells
			const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices2,
			//no sequence indices
	
			//"global" parameters for both WriteToFile instances
			const std::array<bool, 7> filteredValues):
		  
			coupling::AsymmetricalJunctorInterface<dim>( 
				inputCellVector1,
				outputCellVector1,
				mamicoCellIndices1,
				
				inputCellVector2,
				mamicoCellIndices2,
				filteredValues,
		   		"NeuralNet")
		{
			#ifdef DEBUG_NeuralNet
				std::cout<<"Constructed Neural Net"<<std::endl;
			#endif

			using namespace tensorflow;
  			using namespace tensorflow::ops;
	
			auto scope = Scope::NewRootScope();
  			auto x = MatMul(scope, {{1, 1}}, {{41}, {1}});
  			ClientSession session(scope);
  			std::vector<Tensor> outputs;
  			auto status = session.Run({x}, &outputs);
  			TF_CHECK_OK(status);
		
  			std::cout << "Underlying Scalar value -> " << outputs[0].flat<int>()<< std::endl;
			
			
		}
      
		~NeuralNet() {
			
			#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Destroyed NeuralNet instance." << std::endl;
			#endif
		}

    void operator()() {
		#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Now calling operator() on NeuralNet instance." << std::endl;
		#endif
		
		//enters all x velocity values from the outer domain into input_vec
		std::vector<float> input_vec;
		for(unsigned int index = 0; index < coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2.size(); index++) input_vec.push_back((coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2[index]->coupling::datastructures::MacroscopicCell<dim>::getCurrentVelocity)()[0]);
		std::cout<<"\n"<<input_vec.size()<<std::endl;
		
		auto input=VecToTensor(input_vec);
		std::cout<<input.DebugString()<<std::endl;
	}

  private:
	
};

