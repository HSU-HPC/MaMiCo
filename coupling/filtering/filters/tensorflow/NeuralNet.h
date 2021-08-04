/*
Defines the neural network class used for Standalone.cpp
*/
#include "tensorflow/cc/ops/array_ops.h"
#include "tensorflow/cc/framework/ops.h"
#include "tensorflow/cc/ops/state_ops.h"
#include "tensorflow/cc/ops/math_ops.h"

#include <iostream>
#include <map>
#include <fstream>
#include <chrono>
#include <iomanip>
#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/cc/framework/gradients.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/lib/io/path.h"
#include "tensorflow/core/summary/summary_file_writer.h"
#include "tensorflow/cc/tools/freeze_saved_model.h"
#include "tensorflow/core/util/events_writer.h"
#include "tensorflow/cc/framework/scope.h"
#include "TensorCalc.h"

using namespace tensorflow;
using namespace tensorflow::ops;

auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize);

class NeuralNet
{
private:
	//network dimensions
	int input_size, middle_size, output_size;
	
	//scope for network graphs/operations
	Scope net_scope;
    
	//network session pointer
    std::unique_ptr<ClientSession> net_session;
	
	//output objects containing placeholders and network operations
    Output input_placeholder;
    Output label_placeholder;
    Output out_classification;
	Output out_loss_var;
	
    //maps for network weight tensors, their shapes and initializers
    std::map<string, Output> m_vars;
    std::map<string, TensorShape> m_shapes;
    std::map<string, Output> m_assigns;
	
    //vectors used for network initialization/creation
    std::vector<Output> v_weights_biases;
	std::vector<Output> grad_outputs;
	
	//contains adam operations for optimization
    std::vector<Operation> v_out_grads;

public:
	NeuralNet():input_size(0), middle_size(0), output_size(0), net_scope(Scope::NewRootScope()){}
    NeuralNet(int in, int middle, int out): input_size(in), middle_size(middle), output_size(out),  net_scope(Scope::NewRootScope()) {} 
	
	//sets network dimensions
	void setDims(int in, int middle, int out);
	
	//xavier initializer for dense layer weights with given dimensions
    Input Init(Scope scope, int in_chan, int out_chan);
	
	//creates dense layer to add to network graph
	//in_units denote size of previous layer, out_units is the size of this layer
	//bActivation=true adds ReLU activation
	//input is the previous layer object/input placeholder
    Input AddDenseLayer(string idx, Scope scope, int in_units, int out_units, bool bActivation, Input input);
	
	//similar to AddDenseLayer but returns Output object so it can serve as last layer of a graph
    Output AddOutLayer(Scope scope, int in_units, int out_units, Input input);
	
	//creates neural network graph and adds initialization operations to initializer list, this list can be executed with Initialize()
    Status CreateNNGraph();
	
	//creates optimization graph for the network
	//this adds Adam optimizer operations for each layer and also adds the required initialization operations to the initializer list
    Status CreateOptimizationGraph(float learning_rate);
	
	//can be called to change the learning rate of the Adam optimizer even after the network has been initialized
	Status UpdateOptimizationGraph(float learning_rate);
	
	//executes all operations in m_assigns
    Status Initialize();
	
	//trains network graph with provided input and label data, writes prediction to results and mean squared error to loss
    Status Train(Tensor& image_batch, Tensor& label_batch, std::vector<std::vector<float>>& results, float& loss);
	
	//computes network prediction for given input and writes it to results
    Status Predict(Tensor image, std::vector<float>& result);
};

