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

template<unsigned int dim>
class coupling::NeuralNetJunctor : public coupling::AsymmetricalJunctorInterface<dim>{
	using coupling::AsymmetricalJunctorInterface<dim>::_cellIndices2;
	using coupling::FilterInterface<dim>::_outputCells;
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
			const std::array<bool, 7> filteredValues,
			
			std::string nnType,
			int epochs,
			int batchSize,
			float learningRate,
			std::string modelPath):
		  
			coupling::AsymmetricalJunctorInterface<dim>( 
				inputCellVector1,
				outputCellVector1,
				mamicoCellIndices1,
				
				inputCellVector2,
				mamicoCellIndices2,
				filteredValues,
		   		"NeuralNetJunctor"),
			
			nnType_(nnType),
			epochs_(epochs),
			batchSize_(batchSize),
			learningRate_(learningRate),
			modelPath_(modelPath)	
		{
			using namespace tensorflow;
  			using namespace tensorflow::ops;
			
			std::string folder="../../filtering/filters/tensorflow/";		
			if(nnType_=="standalone"){
				
				//loads input and label data from given .csv files into tensors
				//the first two integers in SetupBatchesExcludingGhost denote the columns of the csv that contain the x component of the velocity
				//the last integer is the batch size used for training
				std::vector<tensorflow::Tensor> input, label;
				std::string addendum="_clean";
				SetupBatchesExcludingGhost(input, folder+"clean_data/writer2"+addendum+".csv", 8, label, folder+"clean_data/writer1"+addendum+".csv", 4, batchSize);
				
				//init containers required by NN
				std::vector<std::vector<float>> results;
				std::vector<float> res;
				float loss,total=0.0;
				const int nBatches=input.size();
				
				//init NN
				NN.setDims(input[0].dim_size(1), input[0].dim_size(1), label[0].dim_size(1));
				NN.CreateNNGraph();
				NN.CreateOptimizationGraph(learningRate_);
				NN.Initialize();
				
				
				//train NN with input and label data
				for(int i=0; i<epochs; ++i){
					for(int j=0; j<nBatches; j++){
						NN.Train(input[j], label[j], results, loss); 
						total+=loss;
					}
					std::cout<<"Loss: "<<total/(float)nBatches<<std::endl;
					total=0.0;
					if(i%50==0) std::cout<<"\n"<<"Epoch: "<<i<<"\n"<<std::endl;
				}
			}
			//nnType_ can only be "standalone" or "loadModel", this is checked before constructing NeuralNetJunctor instance
			else{
				//loads SavedModel file into "model"
				SessionOptions session_options = SessionOptions();
				RunOptions run_options = RunOptions();
				Status status = LoadSavedModel(session_options, run_options, modelPath_, {kSavedModelTagServe}, &model);	
				//checks if SavedModel was loaded correctly
				if (!status.ok()) {
					std::cout << "NeuralNetworkJunctor LoadSavedModel Failed: " << status.ToString()<<std::endl;
					exit(2);
				}
				
				//gets names of input and output of the model to later address them
				//names can also be accessed from console with: saved_model_cli show --dir ~/local/tensorflow_cc-master/test/model --all
				inputName_ = model.GetSignatures().at("serving_default").inputs().begin()->second.name();
				outputName_ =model.GetSignatures().at("serving_default").outputs().begin()->second.name();
			}

			#ifdef DEBUG_NeuralNet
				std::cout<<"Constructed Neural Net"<<std::endl;
			#endif
		}
		
		/*
		~NeuralNetJunctor() {
			
			#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Destroyed NeuralNet instance." << std::endl;
			#endif
		}
		*/

    void operator()() {
		#ifdef DEBUG_NeuralNet
				std::cout << "    NeuralNet: Now calling operator() on NeuralNet instance." << std::endl;
		#endif
		
		//enters all x velocity values from the outer domain into inputVec and converts that to inputTensor
		std::vector<float> inputVec;
		for(unsigned int index = 0; index < coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2.size(); index++){ 
			
			if(isNotGhost(_cellIndices2[index])){
				inputVec.push_back((coupling::AsymmetricalJunctorInterface<dim>::_inputCellVector2[index]->coupling::datastructures::MacroscopicCell<dim>::getMicroscopicMomentum)()[0]);
			}
		}
		auto inputTensor=VecToTensor(inputVec);
		
		Tensor resultTensor;
		
		if(nnType_=="standalone"){
			std::vector<float> resultVec;
			
			//generate prediction
			NN.Predict(inputTensor, resultVec);
			resultTensor=VecToTensor(resultVec);
			std::cout<<inputTensor.DebugString()<<std::endl;
			std::cout<<resultTensor.DebugString()<<std::endl;
		}
		else{
				std::vector<Tensor> resultVec;
				
				//generate prediction
				Status runStatus = model.GetSession()->Run({{inputName_, inputTensor}}, {outputName_}, {}, &resultVec);
				if (!runStatus.ok()) {
					std::cout << "NeuralNetworkJunctor Run Failed: " << runStatus.ToString()<<std::endl;
					exit(3);
				}
				resultTensor=resultVec[0];
				std::cout<<inputTensor.DebugString()<<std::endl;
				std::cout<<resultTensor.DebugString()<<std::endl;
		}
		tarch::la::Vector<dim,double> velocity(0);
		
		//the result is written to the currentVelocity vector instead of macroscopicMomentum because calculating the momentum would require macroscopicMass, 
		//which the neural network in its current configuration does not predict
		for(unsigned int index = 0; index < coupling::AsymmetricalJunctorInterface<dim>::_outputCells.size(); index++){ 
			velocity[0]=resultTensor.flat<float>()(index);
			(coupling::AsymmetricalJunctorInterface<dim>::_outputCells[index]->coupling::datastructures::MacroscopicCell<dim>::setCurrentVelocity)(velocity);
		}
	}

  private:
	
	const std::string nnType_;
	int epochs_, batchSize_;
	float learningRate_;
	const std::string modelPath_;
	
	NeuralNet NN;
	SavedModelBundleLite model;
	std::string inputName_, outputName_;
};

//this only works for a 12³ domain
inline bool isNotGhost(tarch::la::Vector<3, unsigned int> indices){
	for(int i=0; i<3;i++){
		if(indices[i]==0 || indices[i]==13) return false;
	}
	return true;
}