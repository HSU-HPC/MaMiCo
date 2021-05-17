
#include "NeuralNet.h"


using namespace std;
using namespace tensorflow;
using namespace tensorflow::ops;

auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize) {
  int32_t* tensor_ptr = tensor.flat<int32_t>().data();
  std::vector<int32_t> v(tensor_ptr, tensor_ptr + tensorSize);
  return v;
}

void NeuralNet::CreateNN(int in, int middle, int out){
	input_size=in;
	middle_size=middle;
	output_size=out;
}

Input NeuralNet::XavierInit(Scope scope, int in_chan, int out_chan)
{
    float std;
    Tensor t;
        std = sqrt(6.f/(in_chan+out_chan));
        Tensor ts(DT_INT64, {2});
        auto v = ts.vec<int64>();
        v(0) = in_chan;
        v(1) = out_chan;
        t = ts;
    auto rand = RandomUniform(scope, t, DT_FLOAT);
    return Multiply(scope, Add(scope, rand, 0.1f), std*2.f);
}

Input NeuralNet::AddDenseLayer(string idx, Scope scope, int in_units, int out_units, bool bActivation, Input input)
{
    TensorShape sp = {in_units, out_units};
    m_vars["W"+idx] = Variable(scope.WithOpName("W"), sp, DT_FLOAT);
    m_shapes["W"+idx] = sp;
    m_assigns["W"+idx+"_assign"] = Assign(scope.WithOpName("W_assign"), m_vars["W"+idx], XavierInit(scope, in_units, out_units));
    sp = {out_units};
    m_vars["B"+idx] = Variable(scope.WithOpName("B"), sp, DT_FLOAT);
    m_shapes["B"+idx] = sp;
    m_assigns["B"+idx+"_assign"] = Assign(scope.WithOpName("B_assign"), m_vars["B"+idx], Input::Initializer(0.f, sp));
    auto dense = Add(scope.WithOpName("Dense_b"), MatMul(scope.WithOpName("Dense_w"), input, m_vars["W"+idx]), m_vars["B"+idx]);
    if(bActivation)
        return Relu(scope.WithOpName("Relu"), dense);
    else
        return dense;
}





/*
Status NeuralNet::FreezeSave(string& file_name)
{
    std::vector<Tensor> out_tensors;
    //Extract: current weights and biases current values
    TF_CHECK_OK(t_session->Run(v_weights_biases , &out_tensors));
    unordered_map<string, Tensor> variable_to_value_map;
    int idx = 0;
    for(Output o: v_weights_biases)
    {
        variable_to_value_map[o.node()->name()] = out_tensors[idx];
        idx++;
    }
    GraphDef graph_def;
    TF_CHECK_OK(t_root.ToGraphDef(&graph_def));
    //call the utility function (modified)
    SavedModelBundle saved_model_bundle;
    SignatureDef signature_def;
    (*signature_def.mutable_inputs())[input_batch_var.name()].set_name(input_batch_var.name());
    (*signature_def.mutable_outputs())[out_classification.name()].set_name(out_classification.name());
    MetaGraphDef* meta_graph_def = &saved_model_bundle.meta_graph_def;
    (*meta_graph_def->mutable_signature_def())["signature_def"] = signature_def;
    *meta_graph_def->mutable_graph_def() = graph_def;
    SessionOptions session_options;
    saved_model_bundle.session.reset(NewSession(session_options));//even though we will not use it
    GraphDef frozen_graph_def;
    std::unordered_set<string> inputs;
    std::unordered_set<string> outputs;
    TF_CHECK_OK(FreezeSavedModel(saved_model_bundle, &frozen_graph_def, &inputs, &outputs));

    //write to file
    return WriteBinaryProto(Env::Default(), file_name, frozen_graph_def);
}
*/

Status NeuralNet::CreateGraphForNN()
{
    input_batch_var = Placeholder(t_root.WithOpName(input_name), DT_FLOAT);

    int in_units = input_size;
    int out_units = middle_size;
    Scope scope_dense1 = t_root.NewSubScope("Dense1_layer");
    auto relu1 = AddDenseLayer("1", scope_dense1, in_units, out_units, true, input_batch_var);
/*
	in_units = out_units;
    out_units = middle_size;
    Scope scope_dense2 = t_root.NewSubScope("Dense2_layer");
    auto relu2 = AddDenseLayer("2", scope_dense2, in_units, out_units, true, relu1);
*/	
    in_units = out_units;
    out_units = output_size;
    Scope scope_dense3 = t_root.NewSubScope("Dense3_layer");
    auto logits = AddDenseLayer("2", scope_dense3, in_units, out_units, false, relu1);

    out_classification =Multiply(t_root.WithOpName(out_name), logits, 1.0f);//Sigmoid(t_root.WithOpName(out_name), logits);
    return t_root.status();
}

Status NeuralNet::CreateOptimizationGraph(float learning_rate)
{
    input_labels_var = Placeholder(t_root.WithOpName("inputL"), DT_FLOAT);
    Scope scope_loss = t_root.NewSubScope("Loss_scope");
    out_loss_var = Mean(scope_loss.WithOpName("Loss"), SquaredDifference(scope_loss, out_classification, input_labels_var), {0});
    TF_CHECK_OK(scope_loss.status());
    for(pair<string, Output> i: m_vars)
        v_weights_biases.push_back(i.second);
    std::vector<Output> grad_outputs;
    TF_CHECK_OK(AddSymbolicGradients(t_root, {out_loss_var}, v_weights_biases, &grad_outputs));
    int index = 0;
    for(pair<string, Output> i: m_vars)
    {
        //Applying Adam
        string s_index = to_string(index);
        auto m_var = Variable(t_root, m_shapes[i.first], DT_FLOAT);
        auto v_var = Variable(t_root, m_shapes[i.first], DT_FLOAT);
        m_assigns["m_assign"+s_index] = Assign(t_root, m_var, Input::Initializer(0.f, m_shapes[i.first]));
        m_assigns["v_assign"+s_index] = Assign(t_root, v_var, Input::Initializer(0.f, m_shapes[i.first]));

        auto adam = ApplyAdam(t_root, i.second, m_var, v_var, 0.f, 0.f, learning_rate, 0.9f, 0.999f, 0.00000001f, {grad_outputs[index]});
        v_out_grads.push_back(adam.operation);
        index++;
    }
    return t_root.status();
}

Status NeuralNet::Initialize()
{
    if(!t_root.ok())
        return t_root.status();
    
    std::vector<Output> ops_to_run;
    for(pair<string, Output> i: m_assigns)
        ops_to_run.push_back(i.second);
    t_session = std::unique_ptr<ClientSession>(new ClientSession(t_root));
    TF_CHECK_OK(t_session->Run(ops_to_run, nullptr));
    /* uncomment if you want visualization of the model graph
    GraphDef graph;
    TF_RETURN_IF_ERROR(t_root.ToGraphDef(&graph));
    SummaryWriterInterface* w;
    TF_CHECK_OK(CreateSummaryFileWriter(1, 0, "/Users/bennyfriedman/Code/TF2example/TF2example/graphs", ".cnn-graph", Env::Default(), &w));
    TF_CHECK_OK(w->WriteGraph(0, make_unique<GraphDef>(graph)));
    */
    return Status::OK();
}

Status NeuralNet::TrainNN(Tensor& image_batch, Tensor& label_batch, std::vector<std::vector<float>>& results, float& loss)
{
    if(!t_root.ok())
        return t_root.status();
    
    std::vector<Tensor> out_tensors;

	
    TF_CHECK_OK(t_session->Run({{input_batch_var, image_batch}, {input_labels_var, label_batch}}, {out_loss_var, out_classification}, v_out_grads, &out_tensors));
	
    loss = tensorMean(out_tensors[0]);
	
	results=TensorToVec(out_tensors[1]);
	/*
    auto mat1 = label_batch.matrix<float>();
    auto mat2 = out_tensors[1].matrix<float>();
    for(int i = 0; i < mat1.dimension(0); i++)
        results.push_back((fabs(mat2(i, 0) - mat1(i, 0)) > 0.5f)? 0 : 1);
	*/
    return Status::OK();
}

Status NeuralNet::Predict(Tensor& image, std::vector<float>& result)
{
    if(!t_root.ok())
        return t_root.status();
    
    std::vector<Tensor> out_tensors;
    //Inputs: image, drop rate 1 and skip drop.
    TF_CHECK_OK(t_session->Run({{input_batch_var, image}}, {out_classification}, &out_tensors));
    result=TensorToVec(out_tensors[0])[0];
    return Status::OK();
}

void write_scalar(tensorflow::EventsWriter* writer, double wall_time, tensorflow::int64 step,
                  const std::string& tag, float simple_value) {
  tensorflow::Event event;
  event.set_wall_time(wall_time);
  event.set_step(step);
  tensorflow::Summary::Value* summ_val = event.mutable_summary()->add_value();
  summ_val->set_tag(tag);
  summ_val->set_simple_value(simple_value);
  writer->WriteEvent(event);
}



