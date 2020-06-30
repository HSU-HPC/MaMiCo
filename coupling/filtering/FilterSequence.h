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
    	FilterSequence(const coupling::IndexConversion<dim>* indexConversion, const char* name):
      	FilterSequence(indexConversion, name, nullptr){}

    	FilterSequence(const coupling::IndexConversion<dim>* indexConversion, const char* name, FilterSequence* input) :
    	_indexConversion(indexConversion), _name(name), _input(input), _isOutput(false){
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Created new sequence named " << _name << ".";
			if(_input) std::cout << " It will use " << _input->getName() << " as input.";
			std::cout << std::endl;
        	#endif
    	}


    	~FilterSequence(){
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Sequence named " << _name << " deconstructed." << std::endl;
        	#endif
    	}

		//interprets specification of filters, initilizes filter objects. Nonzero return implies failure.
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);

    	//Loads cell vector and indices (usually) from input sequence. Not part of constructor because it doesnt get called when sequences are instanced (but in coupling::FilterPipeline<dim>apply(...) instead...)
    	void fillSequenceData(
				const std::vector<coupling::datastructures::MacroscopicCell<dim>* > mdMacroscopicCells,
			   	const std::vector<tarch::la::Vector<dim, unsigned int>> mdCellIndices);

    	const char* getName() { return _name; }
    	const char* getInput() { return _input; }

		bool isOutput() { return _isOutput; };
		void setAsOutput() { _isOutput = true; };

    	void setMacroscopicCells(std::vector<coupling::datastructures::MacroscopicCell<dim>* > macroscopicCells) { _macroscopicCells = macroscopicCells; }
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() { return _macroscopicCells; } //non const
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() const { return _macroscopicCells; } //const 

    	void setCellIndices(const std::vector<tarch::la::Vector<dim, unsigned int>>& cellIndices) { _cellIndices = cellIndices; }
    	std::vector<tarch::la::Vector<dim, unsigned int>> getCellIndices() { return _cellIndices; }
       

		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	
        
	private:
		const coupling::IndexConversion<dim>* _indexConversion;
    	const char* _name;
    	FilterSequence* _input;
		bool _isOutput;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _macroscopicCells;
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;

		std::vector<coupling::FilterInterface<dim> *> _filters;
};

//inlcude implementation
#include "coupling/filtering/FilterSequence.cpph"
