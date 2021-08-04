// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_NeuralNet

#include "coupling/filtering/interfaces/JunctorInterface.h"

#include "tensorflow/NeuralNet.h"
#include "tensorflow/TensorCalc.h"
#include "tensorflow/ConsoleOutput.h"
#include "tensorflow/ReadCSV.h"

#include "tensorflow/cc/ops/array_ops.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/cc/framework/ops.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/cc/ops/state_ops.h"
#include "tensorflow/cc/ops/math_ops.h"
#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"

#include <cstdlib>

bool isNotGhost(tarch::la::Vector<3, unsigned int> indices);

namespace coupling {
  template<unsigned int dim>
  class NeuralNetJunctor;
}

/*auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize) {
  int32_t* tensor_ptr = tensor.flat<int32_t>().data();
  std::vector<int32_t> v(tensor_ptr, tensor_ptr + tensorSize);
  return v;
}
*/

template<unsigned int dim>
class coupling::NeuralNetJunctor : public coupling::AsymmetricalJunctorInterface<dim>{
	public:
		NeuralNetJunctor(
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
		   		"NeuralNetJunctor"),
			inputCellIndices(mamicoCellIndices2)
		{
			using namespace tensorflow;
  			using namespace tensorflow::ops;
			
			//store cell indices to later check for ghost cells
			
			//this switches between the purely TensorFlow C++ and the Python SavedModel approach
			//either standalone or loadModel
			networkType="loadModel";
			
			//this sets the directory containing the required .csv files, needs to be adjusted accordingly
			std::string folder=std::getenv("HOME")+std::string("/local/tensorflow_cc-master/test/");

			//loads config from csv file, contains learning rate and number of epochs
			std::vector<std::string> config=getConfigFromCSV(folder+"config_main.csv");
			
			if(networkType=="standalone"){
				//loads input and label data from given .csv files into tensors
				//the first two integers in SetupBatchesExcludingGhost denote the columns of the csv that contain the x component of the velocity
				//the last integer is the batch size used for training
				std::vector<tensorflow::Tensor> input, label;
				std::string addendum="_clean";
				SetupBatchesExcludingGhost(input, folder+"writer2"+addendum+".csv", 8, label, folder+"writer1"+addendum+".csv", 4, 1);
				
				//init containers required by NN
				std::vector<std::vector<float>> results;
				std::vector<float> res;
				float loss,total=0.0;
				const int nBatches=input.size();
				
				//init NN
				NN.setDims(input[0].dim_size(1), input[0].dim_size(1), label[0].dim_size(1));
				NN.CreateNNGraph();
				NN.CreateOptimizationGraph(std::stof(config[0]));
				NN.Initialize();
				
				
				//train NN with input and label data
				for(int i=0; i<std::stoi(config[1]); ++i){
					for(int j=0; j<nBatches; j++){
						NN.Train(input[j], label[j], results, loss); 
						total+=loss;
					}
					std::cout<<"Loss: "<<total/(float)nBatches<<std::endl;
					total=0.0;
					if(i%100==0) std::cout<<"\n"<<"Epoch: "<<i<<"\n"<<std::endl;
				}
			}
			else if(networkType=="loadModel"){
				//loads SavedModel file into "model"
				SessionOptions session_options = SessionOptions();
				RunOptions run_options = RunOptions();
				Status status = LoadSavedModel(session_options, run_options, folder+config[2], {kSavedModelTagServe}, &model);	
				//checks if SavedModel was loaded correctly
				if (!status.ok()) {
					std::cout << "NeuralNetworkJunctor LoadSavedModel Failed: " << status.ToString()<<std::endl;
					exit(2);
				}
				
				//gets names of input and output of the model to later address them
				//names can also be accessed from console with: saved_model_cli show --dir ~/local/tensorflow_cc-master/test/model --all
				input_name = model.GetSignatures().at("serving_default").inputs().begin()->second.name();
				output_name=model.GetSignatures().at("serving_default").outputs().begin()->second.name();
			}
			else{ 
				std::cout<<"NeuralNetworkJunctor: unknown networkType"<<std::endl;
				exit(1);
			}
			#ifdef DEBUG_NeuralNet
				std::cout<<"Constructed Neural Net"<<std::endl;
			#endif
		}
      
		~NeuralNetJunctor() {
			
			#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Destroyed NeuralNet instance." << std::endl;
			#endif
		}

    void operator()() {
		#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Now calling operator() on NeuralNet instance." << std::endl;
		#endif
		
		//enters all x velocity values from the outer domain into input_vec and converts that to input_pred Tensor
		std::vector<float> input_vec;
		for(unsigned int index = 0; index < coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2.size(); index++){ 
			
			if(isNotGhost(inputCellIndices[index])){
				input_vec.push_back((coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2[index]->coupling::datastructures::MacroscopicCell<dim>::getMicroscopicMomentum)()[0]);
			}
		}
		auto input_tensor=VecToTensor(input_vec);
		
		
		
		if(networkType=="standalone"){
			std::vector<float> result_vec;
			
			//generate prediction
			NN.Predict(input_tensor, result_vec);
			auto result=VecToTensor(result_vec);
			std::cout<<input_tensor.DebugString()<<std::endl;
			std::cout<<result.DebugString()<<std::endl;
		}
		else if(networkType=="loadModel"){
				std::vector<Tensor> result_vec;
				
				//generate prediction
				Status runStatus = model.GetSession()->Run({{input_name, input_tensor}}, {output_name}, {}, &result_vec);
				if (!runStatus.ok()) {
					std::cout << "NeuralNetworkJunctor Run Failed: " << runStatus.ToString()<<std::endl;
					exit(3);
				}
				std::cout<<input_tensor.DebugString()<<std::endl;
				std::cout<<result_vec[0].DebugString()<<std::endl;
		}
		
	}

  private:
	NeuralNet NN;
	SavedModelBundleLite model;
	std::string networkType, input_name, output_name;
	
	//this is a workaround because _cellIndices2 from AsymmetricalJunctorInterface is private
	const std::vector<tarch::la::Vector<dim, unsigned int>> inputCellIndices;
};

//this only works for a 12³ domain
inline bool isNotGhost(tarch::la::Vector<3, unsigned int> indices){
	for(int i=0; i<3;i++){
		if(indices[i]==0 || indices[i]==13) return false;
	}
	return true;
}