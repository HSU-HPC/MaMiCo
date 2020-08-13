// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_PIPELINE
#define POST_MULTI_INSTANCE_FILTERING_YES true
#define POST_MULTI_INSTANCE_FILTERING_NO false
#include "tarch/tinyxml2/tinyxml2.h"
#include "coupling/filtering/FilterSequence.h"
#include "coupling/IndexConversion.h"


/*
 * Manages different branches of filtering sequences.
 * These filtering sequences may be interdependant by using another's sequences input or completely isolated.
 * As this entire filtering process is applied during MD to Macro communication, it uses the MD simulation's output Macro-Cells as input and output.
 * All configuration is made using an XML-config file and does not require recompilation when modified.
 * @Author Felix Maurer
 */
namespace coupling{
    template<unsigned int dim>
    class FilterPipeline;
}


template<unsigned int dim>
class coupling::FilterPipeline{
    public:
        FilterPipeline(
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > mdCells,
			const coupling::IndexConversion<dim>* indexConversion,
			const tarch::utils::MultiMDService<dim>& multiMDService,
			const std::string cfgpath = "filter_pipeline.xml");

        FilterPipeline(
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > mdCells,
			const coupling::IndexConversion<dim>* indexConversion,
			const tarch::utils::MultiMDService<dim>& multiMDService,
			bool postMultiInstance,
			const std::string cfgpath = "filter_pipeline.xml");
               
        ~FilterPipeline() {
            for(auto piSequence : _piSequences) delete piSequence;
            for(auto miSequence : _miSequences) delete miSequence;
            #ifdef DEBUG_FILTER_PIPELINE
            std::cout << "FP: FilterPipeline deconstructed." << std::endl;
            #endif
        }


		/*
		 * Applies each FilterSequence in order of their appearance in the config file.
		 * Ouput of the specified output-FilterSequence will be written to _md2MacroCells.
		 */
        void operator()();



    private:
		/*
		 * Detects errors in XML-config file.
		 */
    	bool configIsValid(tinyxml2::XMLDocument& cfgfile);

		/*
		 * Interprets configuration of sequences and intializes them. Parameters known:
		 *   -"domain-start"/"domain-end": <dim>-Vector (optional, uses entire md2Macro-domain by default)
		 *   -"input": Name of another FilterSequence previously defined (optional, uses MD output (i.e. _md2MacroCells) by default)
		 * Also detects which sequence will be used as output.
		 */
       	int loadSequencesFromXML(tinyxml2::XMLElement* metaNode);

		/*
		 * Chooses a subspace of the cell (and index) input based on what will be transfered to the macro solver.
		 * This subspace is usually called "md2Macro-domain".
		 * */
		int initMd2MacroDomain(std::vector<coupling::datastructures::MacroscopicCell<dim> *> cells);

      
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _md2MacroCells;
		std::vector<tarch::la::Vector<dim, unsigned int>> _md2MacroCellIndices;		

		const coupling::IndexConversion<dim>* _indexConversion;
		const tarch::utils::MultiMDService<dim>& _multiMDService;
	   	bool _postMultiInstance;
      	

		tinyxml2::XMLDocument _config;

		/*
		 * pi = per instance
		 * mi = post multi-instance
		 */
       	std::vector<coupling::FilterSequence<dim> *> _piSequences; 
       	std::vector<coupling::FilterSequence<dim> *> _miSequences;
};


//include implementation of header
#include "FilterPipeline.cpph"
