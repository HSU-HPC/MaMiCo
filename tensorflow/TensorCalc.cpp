

#include "TensorCalc.h"

float tensorMean(tensorflow::Tensor t){
	const int max=t.dim_size(0);
	float sum=0;
	for(int i=0;i<max;++i){
		sum+=t.flat<float>()(i);
	}
	return sum/(float)max;
}

std::vector<float> operator-=(std::vector<float>& vec1, std::vector<float> vec2){
	if(vec1.size()!=vec2.size()){ std::cout<<"\nVectors are not the same size"<<std::endl; return std::vector<float>(1,-1.0f);}
	
	for(int i=0; i<(int)vec1.size();++i){
		vec1[i]-=vec2[i];
	}
	return vec1;
}

float maxInColumn(std::vector<std::vector<float>> vec, int column){
	float max = vec[0][column];
	for (int i=1;i<(int)vec.size();++i){
		if(vec[i][column]>max) max=vec[i][column];
	}
	return max;
}

const tensorflow::Tensor VecToTensor(std::vector<float> vec){
	
	tensorflow::Tensor input(tensorflow::DT_FLOAT, tensorflow::TensorShape({1,(int)vec.size()}));
	auto input_map = input.tensor<float, 2>();
	
	for(int j=0;j<(int)vec.size();++j){
		input_map(0, j) = vec[j];
	}
	return input;
}

std::vector<std::vector<float>> TensorToVec(tensorflow::Tensor t){
	std::vector<std::vector<float>> vec(t.dim_size(0),std::vector<float>(t.dim_size(1), -1));
	auto input_map = t.tensor<float, 2>();
	for(int i=0;i<(int)t.dim_size(0);++i){
		for(int j=0;j<(int)t.dim_size(1);++j){
			 vec[i][j]=input_map(i, j);
		}
	}
	return vec;
}

//returns single 1d tensor out of 2d tensor
tensorflow::Tensor getTensorByIndex(tensorflow::Tensor t, int index){
	return VecToTensor(TensorToVec(t)[index]);
}




