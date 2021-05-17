
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

//using namespace std;
using namespace tensorflow;
using namespace tensorflow::ops;

auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize);

class NeuralNet
{
private:
	int input_size, middle_size, output_size;
    Scope i_root;
    Output file_name_var;
    Output image_tensor_var;
    Scope a_root;
    Output aug_tensor_input;
    Output aug_tensor_output;
    Scope t_root;
    std::unique_ptr<ClientSession> t_session;
    std::unique_ptr<Session> f_session;
    Output input_batch_var;
    string input_name = "input";
    Output input_labels_var;
    Output drop_rate_var; //use real drop rate in training and 1 in validating
    string drop_rate_name = "drop_rate";
    Output skip_drop_var; //use 0 in trainig and 1 in validating
    string skip_drop_name = "skip_drop";
    Output out_classification;
    string out_name = "output_classes";
    Output logits;
    //Network maps
    std::map<string, Output> m_vars;
    std::map<string, TensorShape> m_shapes;
    std::map<string, Output> m_assigns;
    //Loss variables
    std::vector<Output> v_weights_biases;
    std::vector<Operation> v_out_grads;
    Output out_loss_var;
    InputList MakeTransforms(int batch_size, Input a0, Input a1, Input a2, Input b0, Input b1, Input b2);
public:
	NeuralNet():input_size(0), middle_size(0), output_size(0), i_root(Scope::NewRootScope()), a_root(Scope::NewRootScope()), t_root(Scope::NewRootScope()){}
    NeuralNet(int in, int middle, int out): input_size(in), middle_size(middle), output_size(out), i_root(Scope::NewRootScope()), a_root(Scope::NewRootScope()), t_root(Scope::NewRootScope()) {} 
	void CreateNN(int in, int middle, int out);
    Input XavierInit(Scope scope, int in_chan, int out_chan);
    Input AddDenseLayer(string idx, Scope scope, int in_units, int out_units, bool bActivation, Input input);
    Status CreateGraphForNN();
    Status CreateOptimizationGraph(float learning_rate);
    Status Initialize();
    Status TrainNN(Tensor& image_batch, Tensor& label_batch, std::vector<std::vector<float>>& results, float& loss);
    Status ValidateNN(Tensor& image_batch, Tensor& label_batch, std::vector<float>& results);
    Status Predict(Tensor& image, std::vector<float>& result);
    Status FreezeSave(string& file_name);
    Status LoadSavedModel(string& file_name);
    Status PredictFromFrozen(Tensor& image, int& result);
    Status CreateAugmentGraph(int batch_size, int image_side, float flip_chances, float max_angles, float sscale_shift_factor);
    Status RandomAugmentBatch(Tensor& image_batch, Tensor& augmented_batch);
    Status WriteBatchToImageFiles(Tensor& image_batch, string folder_name, string image_name);
};


void write_scalar(tensorflow::EventsWriter* writer, double wall_time, tensorflow::int64 step, const std::string& tag, float simple_value);

//#include "NN.cpph"