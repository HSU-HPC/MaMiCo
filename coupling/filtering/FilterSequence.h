// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include "coupling/IndexConversion.h"
#include "tarch/configuration/ParseConfiguration.h"

//INCLUDE ALL FILTERS HERE
#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/filters/Gauss.h"

//Filter Sequences are used to logically group filters that will be used in chronological order.
//@Author Felix Maurer

namespace coupling{
    template<unsigned int dim>
    class FilterSequence;
}

template<unsigned int dim>
class coupling::FilterSequence {
	public:
    	FilterSequence(const coupling::IndexConversion<dim>* indexConversion, const char* name, std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices) :
    	_indexConversion(indexConversion), _name(name), _cellIndices(cellIndices), _isOutput(false){
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Created new sequence named " << _name << "." << std::endl;
        	#endif
    	}


    	~FilterSequence(){
			for (auto v1 : _cellVector1) delete v1;
			for (auto v2 : _cellVector2) delete v2;
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Sequence named " << _name << " deconstructed." << std::endl;
        	#endif
    	}


		//TODO: move to constructor (?)
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);

    	void fillSequenceData(const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& inputCellVector);
		
		//TODO: besprechen
		void updateCellVectors(){
			for(unsigned int index = 0; index < _inputCellVector.size(); index++){
				*(_cellVector1[index]) = *(_inputCellVector[index]);
				*(_cellVector2[index]) = *(_inputCellVector[index]);
			}
		}

    	const char* getName() { return _name; }

		bool isOutput() { return _isOutput; };
		void setAsOutput() { _isOutput = true; };

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() { return _outputCellVector; } //non const
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() const { return _outputCellVector; } //const

    	void setCellIndices(const std::vector<tarch::la::Vector<dim, unsigned int>>& cellIndices) { _cellIndices = cellIndices; }
    	std::vector<tarch::la::Vector<dim, unsigned int>> getCellIndices() { return _cellIndices; }
       
		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	
        
	private:
		const coupling::IndexConversion<dim>* _indexConversion;
    	const char* _name;

		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector;//pointers to macro cells of this sequence's input
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1;//allocated for this sequence only
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector2;//allocated for this sequence only
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCellVector;//pointers to either _cellVector1 or _cellVector2. Use this as output if _isOutput == true and as input for other sequences
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;//all of the above use the same indexing

		bool _isOutput;
		std::vector<coupling::FilterInterface<dim> *> _filters;
};

//inlcude implementation
#include "coupling/filtering/FilterSequence.cpph"
