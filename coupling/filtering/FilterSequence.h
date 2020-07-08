// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include "coupling/IndexConversion.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "coupling/filtering/WriteToFile.h"

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
			for (auto i : _inputCellVector) delete i;
			for (auto o : _outputCellVector) delete o;
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Sequence named " << _name << " deconstructed." << std::endl;
        	#endif
    	}


		//TODO: move to constructor, "md" should be "input" (which in some cases is the same ofc)
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);
    	void fillSequenceData(const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& mdMacroscopicCells);


    	const char* getName() { return _name; }
    	const char* getInput() { return _input; }

		bool isOutput() { return _isOutput; };
		void setAsOutput() { _isOutput = true; };

    	//void setMacroscopicCells(std::vector<coupling::datastructures::MacroscopicCell<dim>* > macroscopicCells) { _macroscopicCells = macroscopicCells; }
		
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() { return _outputCellVector; } //non const
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() const { return _outputCellVector; } //const

    	void setCellIndices(const std::vector<tarch::la::Vector<dim, unsigned int>>& cellIndices) { _cellIndices = cellIndices; }
    	std::vector<tarch::la::Vector<dim, unsigned int>> getCellIndices() { return _cellIndices; }
       
		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	
        
	private:
		const coupling::IndexConversion<dim>* _indexConversion;
    	const char* _name;
    	FilterSequence* _input;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector;
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCellVector;
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;

		bool _isOutput;
		std::vector<coupling::FilterInterface<dim> *> _filters;
};

//inlcude implementation
#include "coupling/filtering/FilterSequence.cpph"
