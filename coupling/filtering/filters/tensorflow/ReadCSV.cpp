#include "ReadCSV.h"

std::string_view CSVRow::operator[](std::size_t index) const{
    return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
}

std::size_t CSVRow::size() const{
    return m_data.size() - 1;
}

void CSVRow::readNextRow(std::istream& str){
	std::getline(str, m_line);

    m_data.clear();
    m_data.emplace_back(-1);
    std::string::size_type pos = 0;
    while((pos = m_line.find(';', pos)) != std::string::npos){
        m_data.emplace_back(pos);
        ++pos;
    }
    // This checks for a trailing semicolon with no data after it.
    pos   = m_line.size();
    m_data.emplace_back(pos);
}


std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   



CSVIterator::CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }

CSVIterator::CSVIterator()                   :m_str(NULL) {}


CSVIterator& CSVIterator::operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}

CSVRow const& CSVIterator::operator*()   const       {return m_row;}

CSVRow const* CSVIterator::operator->()  const       {return &m_row;}

bool CSVIterator::operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}

bool CSVIterator::operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}

std::vector<std::string> getConfigFromCSV(std::string filename){
	std::ifstream file(filename);
	CSVIterator loop(file);
	std::vector<std::string> vec;
	std::string temp;

	for(int i=0;i<(int)(*loop).size();++i){
		temp=std::string((*loop)[i]);
		temp.erase(remove_if(temp.begin(), temp.end(), isspace), temp.end());
		if((*loop)[i]!=""&&(*loop)[i]!=" ")vec.push_back(temp);
	}
	return vec;
}

std::vector<std::vector<float>> getCSVasVec(std::string filename){
	std::ifstream file(filename);
	CSVIterator loop(file);
	std::vector<std::vector<float>> vec;
	std::vector<float> temp;

    for(; loop != CSVIterator(); ++loop){
		temp.clear();
		//this avoids empty lines
		if((*loop)[0]==""||(*loop)[0]==" ") continue;
		for(int i=0;i<(int)(*loop).size();++i){
			//adds all non-empty elements
			if((*loop)[i]!=""&&(*loop)[i]!=" ")temp.push_back((float)std::stod(std::string((*loop)[i])));
		}
		vec.push_back(temp);
	}
	return vec;
}

std::vector<std::vector<float>> getCSVasVec(std::vector<std::string> filenameVec){
	std::vector<std::vector<float>> vec, temp;

    	for(auto filename : filenameVec){
		temp=getCSVasVec(filename);
		vec.insert(vec.end(), temp.begin(), temp.end());
	}
	
	return vec;
}

std::vector<std::vector<float>> getCSVasVec(std::string filename, int index){
	std::vector<std::vector<float>> temp=getCSVasVec(filename);
	std::vector<std::vector<float>> vec(maxInColumn(temp, 0),std::vector<float>());
	for(auto element : temp) vec[element[0]-1].push_back(element[index]);
	return vec;
}

std::vector<std::vector<float>> getCSVasVec(std::string filename, int index, int divisor){
	std::vector<std::vector<float>> temp=getCSVasVec(filename);
	std::vector<std::vector<float>> vec(maxInColumn(temp, 0),std::vector<float>());
	for(auto element : temp) vec[element[0]-1].push_back(element[index]/element[divisor]);
	return vec;
}

std::vector<std::vector<float>> getCSVasVecExcludingGhost(std::string filename, int index){
	std::vector<std::vector<float>> temp=getCSVasVec(filename);
	std::vector<std::vector<float>> vec(maxInColumn(temp, 0),std::vector<float>());
	for(auto element : temp){
		if(element[1]!=0 && element[1]!=13 && element[2]!=0 && element[2]!=13 && element[3]!=0 && element[3]!=13)vec[element[0]-1].push_back(element[index]);
	}
	return vec;
}

std::vector<std::vector<float>> getCSVasVecExcludingGhost(std::string filename, int index, int divisor){
	std::vector<std::vector<float>> temp=getCSVasVec(filename);
	std::vector<std::vector<float>> vec(maxInColumn(temp, 0),std::vector<float>());
	for(auto element : temp){
		if(element[1]!=0 && element[1]!=13 && element[2]!=0 && element[2]!=13 && element[3]!=0 && element[3]!=13)vec[element[0]-1].push_back(element[index]/element[divisor]);
	}
	return vec;
}

const std::vector<tensorflow::Tensor> getCSVasVecOfBatches(std::string filename, int index, int batch_size){
	std::vector<tensorflow::Tensor> vec;
	std::vector<std::vector<float>> temp=getCSVasVec(filename, index);
	std::vector<std::vector<float>> temp2;
	int i=0;
	
	for(std::vector<float> entry : temp){
		i++;
		temp2.push_back(entry);
		if(i%batch_size==0){
			vec.push_back(VecToTensor(temp2));
			temp2.clear();
		}
	}
	
	return vec;
}

const std::vector<tensorflow::Tensor> getCSVasVecOfBatches(std::string filename, int index, int batch_size, int divisor){
	std::vector<tensorflow::Tensor> vec;
	std::vector<std::vector<float>> temp=getCSVasVec(filename, index, divisor);
	std::vector<std::vector<float>> temp2;
	int i=0;
	
	for(std::vector<float> entry : temp){
		i++;
		temp2.push_back(entry);
		if(i%batch_size==0){
			vec.push_back(VecToTensor(temp2));
			temp2.clear();
		}
	}
	
	return vec;
}

const std::vector<tensorflow::Tensor> getCSVasVecOfBatchesExcludingGhost(std::string filename, int index, int batch_size){
	std::vector<tensorflow::Tensor> vec;
	std::vector<std::vector<float>> temp=getCSVasVecExcludingGhost(filename, index);
	std::vector<std::vector<float>> temp2;
	int i=0;
	
	for(std::vector<float> entry : temp){
		i++;
		temp2.push_back(entry);
		if(i%batch_size==0){
			vec.push_back(VecToTensor(temp2));
			temp2.clear();
		}
	}
	
	return vec;
}

const std::vector<tensorflow::Tensor> getCSVasVecOfBatchesExcludingGhost(std::string filename, int index, int batch_size, int divisor){
	std::vector<tensorflow::Tensor> vec;
	std::vector<std::vector<float>> temp=getCSVasVecExcludingGhost(filename, index, divisor);
	std::vector<std::vector<float>> temp2;
	int i=0;
	
	for(std::vector<float> entry : temp){
		i++;
		temp2.push_back(entry);
		if(i%batch_size==0){
			vec.push_back(VecToTensor(temp2));
			temp2.clear();
		}
	}
	
	return vec;
}

void VecToCSV(std::vector<std::vector<float>> vec2d){
	if(std::remove("output.csv")==0) std::cout<<"output.csv was deleted"<<std::endl;
	std::ofstream myfile("output.csv");
    
	for(std::vector<float> vec : vec2d){
		for(float number : vec){
			myfile << number << '\t';
			myfile << "," ;
		}
		myfile << std::endl;
	}
}

void VecToCSV(std::vector<float> vec){
	if(std::remove("output.csv")==0) std::cout<<"output.csv was deleted"<<std::endl;
	std::ofstream myfile("output.csv");
    
    for(float number : vec){
		myfile << number << '\t';
		myfile << "," ;
    }
}


void SetupBatches(std::vector<tensorflow::Tensor>& in, std::string in_path, int in_column, std::vector<tensorflow::Tensor>& label, std::string label_path, int label_column, int batch_size){
	in=getCSVasVecOfBatches(in_path, in_column, batch_size);
	label=getCSVasVecOfBatches(label_path, label_column, batch_size);
}

void SetupBatches(std::vector<tensorflow::Tensor>& in, std::string in_path, int in_column, std::vector<tensorflow::Tensor>& label, std::string label_path, int label_column, int batch_size, int divisor){
	in=getCSVasVecOfBatches(in_path, in_column, batch_size);
	label=getCSVasVecOfBatches(label_path, label_column, batch_size, divisor);
}

void SetupBatchesExcludingGhost(std::vector<tensorflow::Tensor>& in, std::string in_path, int in_column, std::vector<tensorflow::Tensor>& label, std::string label_path, int label_column, int batch_size){
	in=getCSVasVecOfBatchesExcludingGhost(in_path, in_column, batch_size);
	label=getCSVasVecOfBatches(label_path, label_column, batch_size);
}

void SetupBatchesExcludingGhost(std::vector<tensorflow::Tensor>& in, std::string in_path, int in_column, std::vector<tensorflow::Tensor>& label, std::string label_path, int label_column, int batch_size, int divisor){
	in=getCSVasVecOfBatchesExcludingGhost(in_path, in_column, batch_size);
	label=getCSVasVecOfBatches(label_path, label_column, batch_size, divisor);
}