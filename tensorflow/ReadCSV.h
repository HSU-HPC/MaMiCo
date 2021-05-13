/*
see https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c?page=1&tab=votes#tab-top
*/
#pragma once

#include "tensorflow/core/framework/tensor.h"
#include "ConsoleOutput.h"
#include "TensorCalc.h"
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <string_view>



std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str);

class CSVRow
{
    public:
        std::string_view operator[](std::size_t index) const;
        std::size_t size() const;
        void readNextRow(std::istream& str);
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data);

class CSVIterator
{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str);
        CSVIterator();

        // Pre Increment
        CSVIterator& operator++();
        // Post increment
        CSVIterator operator++(int);
        CSVRow const& operator*()   const; 
        CSVRow const* operator->()  const;

        bool operator==(CSVIterator const& rhs);
        bool operator!=(CSVIterator const& rhs);
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

//turns 2d vector into 2d tensor with the same entries
//all 1d vectors must be the same length 
const tensorflow::Tensor VecToTensor(std::vector<std::vector<float>> vec);

std::vector<std::vector<float>> getCSVasVec(std::string filename);

std::vector<std::vector<float>> getCSVasVec(std::vector<std::string> filenameVec);

//returns 2d vector containing the [index]-entry of every line of the csv
//first vector index denotes the timestep-1, second index denotes macro coordinate(from 4,4,4 over 5,4,4 to 10,10,10)
std::vector<std::vector<float>> getCSVasVec(std::string filename, int index);

std::vector<float> getCSVEntryasVec(std::string filename, int index);

std::vector<std::vector<float>> getCSVEntryasVec(std::vector<std::string> filenameVec, int index);

const tensorflow::Tensor getCSVasTensor(std::string filename);

const tensorflow::Tensor getCSVasTensor(std::vector<std::string> filenameVec);

const tensorflow::Tensor getCSVEntriesasTensor(std::vector<std::string> filenameVec, int index);

const tensorflow::Tensor getCSVasTensor(std::string filename, int index);

//#include"ReadCSV.cpph"







