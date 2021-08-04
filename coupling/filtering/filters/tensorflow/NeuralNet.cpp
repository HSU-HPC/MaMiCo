#include "NeuralNet.h"


using namespace std;
using namespace tensorflow;
using namespace tensorflow::ops;

auto tf_tensor_to_vector(tensorflow::Tensor tensor, int32_t tensorSize) {
  int32_t* tensor_ptr = tensor.flat<int32_t>().data();
  std::vector<int32_t> v(tensor_ptr, tensor_ptr + tensorSize);
  return v;
}

void NeuralNet::setDims(int in, int middle, int out){
	input_size=in;
	middle_size=middle;
	output_size=out;
}

Input NeuralNet::Init(Scope scope, int in_chan, int out_chan)
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
    m_assigns["W"+idx+"_assign"] = Assign(scope.WithOpName("W_assign"), m_vars["W"+idx], Init(scope, in_units, out_units));
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

Output NeuralNet::AddOutLayer(Scope scope, int in_units, int out_units, Input input)
{
	std::string idx="2";
    TensorShape sp = {in_units, out_units};
    m_vars["W"+idx] = Variable(scope.WithOpName("W"), sp, DT_FLOAT);
    m_shapes["W"+idx] = sp;
    m_assigns["W"+idx+"_assign"] = Assign(scope.WithOpName("W_assign"), m_vars["W"+idx], Init(scope, in_units, out_units));
    sp = {out_units};
    m_vars["B"+idx] = Variable(scope.WithOpName("B"), sp, DT_FLOAT);
    m_shapes["B"+idx] = sp;
    m_assigns["B"+idx+"_assign"] = Assign(scope.WithOpName("B_assign"), m_vars["B"+idx], Input::Initializer(0.f, sp));
    auto out = Add(scope.WithOpName("Dense_b"), MatMul(scope.WithOpName("Dense_w"), input, m_vars["W"+idx]), m_vars["B"+idx]);
    return out;
}


Status NeuralNet::CreateNNGraph()
{

    input_placeholder = Placeholder(net_scope.WithOpName("input"), DT_FLOAT);
	
    int in_units = input_size;
    int out_units = middle_size;
    Scope scope_dense1 = net_scope.NewSubScope("Dense1_layer");
    auto relu = AddDenseLayer("1", scope_dense1, in_units, out_units, true, input_placeholder);

	in_units = out_units;
    out_units = output_size;
    Scope scope_dense2 = net_scope.NewSubScope("Dense2_layer");
	out_classification=AddOutLayer(scope_dense2, in_units, out_units, relu);

    return net_scope.status();

}

Status NeuralNet::CreateOptimizationGraph(float learning_rate)
{
    label_placeholder = Placeholder(net_scope.WithOpName("inputL"), DT_FLOAT);
    Scope scope_loss = net_scope.NewSubScope("Loss_scope");
    out_loss_var = Mean(scope_loss.WithOpName("Loss"), SquaredDifference(scope_loss, out_classification, label_placeholder), {0});
    TF_CHECK_OK(scope_loss.status());
    for(pair<string, Output> i: m_vars)
        v_weights_biases.push_back(i.second);
    std::vector<Output> grad_outputs;
    TF_CHECK_OK(AddSymbolicGradients(net_scope, {out_loss_var}, v_weights_biases, &grad_outputs));
    int index = 0;
    for(pair<string, Output> i: m_vars)
    {
        string s_index = to_string(index);
        auto m_var = Variable(net_scope, m_shapes[i.first], DT_FLOAT);
        auto v_var = Variable(net_scope, m_shapes[i.first], DT_FLOAT);
        m_assigns["m_assign"+s_index] = Assign(net_scope, m_var, Input::Initializer(0.f, m_shapes[i.first]));
        m_assigns["v_assign"+s_index] = Assign(net_scope, v_var, Input::Initializer(0.f, m_shapes[i.first]));

        auto adam = ApplyAdam(net_scope, i.second, m_var, v_var, 0.f, 0.f, learning_rate, 0.9f, 0.999f, 0.00000001f, {grad_outputs[index]});
        v_out_grads.push_back(adam.operation);
        index++;
    }
    return net_scope.status();
}

Status NeuralNet::UpdateOptimizationGraph(float learning_rate)
{

	v_out_grads.clear();
	std::map<string, Output> m_assigns_new;
	int index = 0;
    for(pair<string, Output> i: m_vars)
    {
        string s_index = to_string(index);
        auto m_var = Variable(net_scope, m_shapes[i.first], DT_FLOAT);
        auto v_var = Variable(net_scope, m_shapes[i.first], DT_FLOAT);
        m_assigns_new["m_assign"+s_index] = Assign(net_scope, m_var, Input::Initializer(0.f, m_shapes[i.first]));
        m_assigns_new["v_assign"+s_index] = Assign(net_scope, v_var, Input::Initializer(0.f, m_shapes[i.first]));

        auto adam = ApplyAdam(net_scope, i.second, m_var, v_var, 0.f, 0.f, learning_rate, 0.9f, 0.999f, 0.00000001f, {grad_outputs[index]});
        v_out_grads.push_back(adam.operation);
        index++;
    }
	
	std::vector<Output> ops_to_run;
    for(pair<string, Output> i: m_assigns_new) ops_to_run.push_back(i.second);

    TF_CHECK_OK(net_session->Run(ops_to_run, nullptr));
	
    return net_scope.status();
}

Status NeuralNet::Initialize()
{
    if(!net_scope.ok())
        return net_scope.status();
    
    std::vector<Output> ops_to_run;
    for(pair<string, Output> i: m_assigns)
        ops_to_run.push_back(i.second);
    net_session = std::unique_ptr<ClientSession>(new ClientSession(net_scope));
    TF_CHECK_OK(net_session->Run(ops_to_run, nullptr));

    return Status::OK();
}

Status NeuralNet::Train(Tensor& image_batch, Tensor& label_batch, std::vector<std::vector<float>>& results, float& loss)
{
    if(!net_scope.ok())
        return net_scope.status();
    
    std::vector<Tensor> out_tensors;

    TF_CHECK_OK(net_session->Run({{input_placeholder, image_batch}, {label_placeholder, label_batch}}, {out_loss_var, out_classification}, v_out_grads, &out_tensors));
	
    loss = tensorMean(out_tensors[0]);
	results=TensorToVec(out_tensors[1]);

    return Status::OK();
}

Status NeuralNet::Predict(Tensor image, std::vector<float>& result)
{
    if(!net_scope.ok())
        return net_scope.status();
    
    std::vector<Tensor> out_tensors;
    TF_CHECK_OK(net_session->Run({{input_placeholder, image}}, {out_classification}, &out_tensors));
    result=TensorToVec(out_tensors[0])[0];
    return Status::OK();
}




