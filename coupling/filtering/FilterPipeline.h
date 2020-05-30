// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <iostream>
#include "tarch/tinyxml2/tinyxml2.h"
#include "coupling/filtering/WriteToFile.h"

//#define DEBUG_FILTER_PIPELINE

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

template<unsigned int dim, class CellService>//template für LinkedCell?
class coupling::FilterSequence {
public:
    FilterSequence(const char* name):
        _name(name), _input(nullptr)

    {
        #ifdef DEBUG_FILTER_PIPELINE
        std::cout << "FP: Created new sequence named " << _name << ". It will use default (md) input." << std::endl;
        #endif
    }
    FilterSequence(const char* name, FilterSequence* input /*, const std::vector<coupling::datastructures::MacroscopicCell<dim>* > macroscopicCells, const unsigned int * const cellIndices*/) :
    _name(name), _input(input)/*, _macroscopicCells(macroscopicCells), _cellIndices(cellIndices)*/
    {
        #ifdef DEBUG_FILTER_PIPELINE
        std::cout << "FP: Created new sequence named " << _name << ". It will use " << _input->getName() << " as input." << std::endl;
        #endif
    }


    ~FilterSequence(){
        #ifdef DEBUG_FILTER_PIPELINE
        std::cout << "FP: Sequence named " << _name << " deconstructed." << std::endl;
        #endif
    }

    //Loads cell vector and indices from input sequence. Not part of constructor because it doesnt get called when sequences are instanced (but in coupling::FilterPipeline<dim>apply(...) instead...)
    void initSequence(std::vector<coupling::datastructures::MacroscopicCell<dim>* > mdMacroscopicCells, std::vector<int> mdCellIndices);

    const char* getName() { return _name; }
    const char* getInput() { return _input; }


    void setMacroscopicCells(std::vector<coupling::datastructures::MacroscopicCell<dim>* > macroscopicCells) { _macroscopicCells = macroscopicCells; }
    std::vector<coupling::datastructures::MacroscopicCell<dim>* > getMacroscopicCells() { return _macroscopicCells; } 
        
    void setCellIndices(std::vector<unsigned int> cellIndices) { _cellIndices = cellIndices; }
    std::vector<unsigned int> getCellIndices() { return _cellIndices; }
        
        
private:
    const char* _name;
    FilterSequence* _input;

    std::vector<coupling::datastructures::MacroscopicCell<dim>> _macroscopicCells;
    std::vector<unsigned int> _cellIndices;

};



template<unsigned int dim, class CellService>
class coupling::FilterPipeline{
    public:
        FilterPipeline(const std::string cfgpath = "filter_pipeline.xml") {
            FilterPipeline(POST_MULTI_INSTANCE_FILTERING_NO, cfgpath);
        }

        FilterPipeline(bool postMultiInstance, const std::string cfgpath = "filter_pipeline.xml"){
            //check if provided file is written in proper XML
            if(_config.LoadFile(cfgpath.c_str()) != tinyxml2::XML_NO_ERROR){
                std::cout << "ERROR: Could not read config for Filter-Pipeline: XML syntax error." << std::endl;
                exit(EXIT_FAILURE);
            }
            
            //check for structural errors in config file
            if (!configIsValid(_config)) exit(EXIT_FAILURE);

            //load sequences
            if (piInitializeSequences(_config.FirstChildElement("filter-pipeline")->FirstChildElement("per-instance"))) exit(EXIT_FAILURE);
            if (postMultiInstance){
                if(miInitializeSequences(_config.FirstChildElement("filter-pipeline")->FirstChildElement("post-multi-instance"))) exit(EXIT_FAILURE);
            }
        }
        
        ~FilterPipeline() {
            for(auto piSequence : _piSequences) delete piSequence;
            for(auto miSequence : _miSequences) delete miSequence;
            #ifdef DEBUG_FILTER_PIPELINE
            std::cout << "FP: FilterPipeline deconstructed." << std::endl;
            #endif
        }

        //apply() is overloaded to distinct being called from a single-instance or multimd purpose

        void apply(CellService* cellService);

        



    private:
       bool configIsValid(tinyxml2::XMLDocument& cfgfile);
       int piInitializeSequences(tinyxml2::XMLElement* perInstanceNode);
       int miInitializeSequences(tinyxml2::XMLElement* multiInstanceNode);
       
       tinyxml2::XMLDocument _config;

       std::vector<coupling::FilterSequence<dim, CellService> *> _piSequences; 
       std::vector<coupling::FilterSequence<dim, CellService> *> _miSequences;
};


//include implementation of header
#include "FilterPipeline.cpph"


#ifdef DEBUG_FILTER_PIPELINE
//this is supposed to be a header file. for debugging purposes only.
int main (void) {
    coupling::FilterPipeline<3> fp = coupling::FilterPipeline<3>(POST_MULTI_INSTANCE_FILTERING_NO);
    fp.apply();
    return 0;
}
#endif


/**
 * TODO:
 *
 * */
