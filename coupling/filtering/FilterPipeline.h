// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include "tarch/tinyxml2/tinyxml2.h"
#include "coupling/filtering/FilterSequence.h"

#define DEBUG_FILTER_PIPELINE

#define POST_MULTI_INSTANCE_FILTERING_YES true
#define POST_MULTI_INSTANCE_FILTERING_NO false

//A pipeline-like approach used for chanining and branching of filters for MacrosCopicCells.
//In genreal "pi" stands for per-istance filtering and "mi" for (post-)multi-instance filtering.
//@Author Felix Maurer

namespace coupling{
    template<unsigned int dim>
    class FilterPipeline;
}


template<unsigned int dim>
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


		template<class CellService>
		//globalVectorCellIndices only means cells possibly relevant for filtering, that is cells in MD domain
        void apply(CellService* cellService);



    private:
    	bool configIsValid(tinyxml2::XMLDocument& cfgfile);
       	int loadSequencesFromXML(tinyxml2::XMLElement* metaNode);
       
      	tinyxml2::XMLDocument _config;

	   	bool _postMultiInstance;

       	std::vector<coupling::FilterSequence<dim> *> _piSequences; 
       	std::vector<coupling::FilterSequence<dim> *> _miSequences;
};


//include implementation of header
#include "FilterPipeline.cpph"



/**
 * TODO:
 *
 * */
