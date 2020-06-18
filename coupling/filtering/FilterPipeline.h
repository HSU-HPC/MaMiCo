// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <iostream>
#include "tarch/tinyxml2/tinyxml2.h"
#include "coupling/filtering/WriteToFile.h"
#include "coupling/filtering/FilterInterface.h"
#include "coupling/services/MacroscopicCellService.h"

#define DEBUG_FILTER_PIPELINE

#define POST_MULTI_INSTANCE_FILTERING_YES true
#define POST_MULTI_INSTANCE_FILTERING_NO false

//A pipeline-like approach used for chanining and branching of filters for MacrosCopicCells.
//In genreal "pi" stands for per-istance filtering and "mi" for (post-)multi-instance filtering.
//@Author Felix Maurer

namespace coupling{
    template<unsigned int dim, class CellService>
    class FilterPipeline;

    template<unsigned int dim, class CellService>
    class FilterSequence;
}

template<unsigned int dim, class CellService>
class coupling::FilterSequence {
	public:
    	FilterSequence(const char* name):
      	FilterSequence(name, nullptr){}

    	FilterSequence(const char* name, FilterSequence* input) :
    	_name(name), _input(input), _isOutput(false){
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Created new sequence named " << _name << ".";
			if(_input) std::cout << " It will use " << _input->getName() << " as input.";
			std::cout << std::endl;
        	#endif
    	}


    	~FilterSequence(){
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FP: Sequence named " << _name << " deconstructed." << std::endl;
        	#endif
    	}

		//interprets specification of filters, initilizes filter objects. Nonzero return implies failure.
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);

    	//Loads cell vector and indices (usually) from input sequence. Not part of constructor because it doesnt get called when sequences are instanced (but in coupling::FilterPipeline<dim>apply(...) instead...)
    	void fillSequenceData(std::vector<coupling::datastructures::MacroscopicCell<dim>* > mdMacroscopicCells, std::vector<int> mdCellIndices);

    	const char* getName() { return _name; }
    	const char* getInput() { return _input; }

		bool isOutput() { return _isOutput; };
		void setAsOutput() { _isOutput = true; };

    	void setMacroscopicCells(std::vector<coupling::datastructures::MacroscopicCell<dim>* > macroscopicCells) { _macroscopicCells = macroscopicCells; }
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() { return _macroscopicCells; } //non const
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getMacroscopicCells() const { return _macroscopicCells; } //const 

    	void setCellIndices(std::vector<unsigned int> cellIndices) { _cellIndices = cellIndices; }
    	std::vector<unsigned int> getCellIndices() { return _cellIndices; }
       

		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	
        
	private:
    	const char* _name;
    	FilterSequence* _input;
		bool _isOutput;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _macroscopicCells;
    	std::vector<unsigned int> _cellIndices;

		std::vector<coupling::FilterInterface<dim> *> _filters;
};



template<unsigned int dim, class CellService>
class coupling::FilterPipeline{
    public:
        FilterPipeline(const std::string cfgpath = "filter_pipeline.xml");

        FilterPipeline(
			bool postMultiInstance,
			const std::string cfgpath = "filter_pipeline.xml");
               
        ~FilterPipeline() {
            for(auto piSequence : _piSequences) delete piSequence;
            for(auto miSequence : _miSequences) delete miSequence;
            #ifdef DEBUG_FILTER_PIPELINE
            std::cout << "FP: FilterPipeline deconstructed." << std::endl;
            #endif
        }


        void apply(const unsigned int * const globalCellIndices, CellService* cellService);



    private:
    	bool configIsValid(tinyxml2::XMLDocument& cfgfile);
       	int loadSequencesFromXML(tinyxml2::XMLElement* metaNode);
       
      	tinyxml2::XMLDocument _config;

	   	bool _postMultiInstance;

       	std::vector<coupling::FilterSequence<dim, CellService> *> _piSequences; 
       	std::vector<coupling::FilterSequence<dim, CellService> *> _miSequences;
};


//include implementation of header
#include "FilterPipeline.cpph"



/**
 * TODO:
 *
 * */
