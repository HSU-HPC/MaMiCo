// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

//#define DEBUG_FILTER_PIPELINE
#define POST_MULTI_INSTANCE_FILTERING_YES true
#define POST_MULTI_INSTANCE_FILTERING_NO false

#include "tarch/tinyxml2/tinyxml2.h"
#include "coupling/filtering/sequencing/FilterSequence.h"
#include "coupling/filtering/sequencing/FilterJunction.h"
#include "coupling/filtering/sequencing/AsymmetricalFilterJunction.h"
#include "coupling/IndexConversionMD2Macro.h"

/*
 * TODO: rework comment
 *
 * Manages different branches of filtering sequences.
 * These filtering sequences may be interdependant by using another's sequences input or completely isolated.
 * As this entire filtering process is applied during MD to Macro communication, it uses the MD simulation's output Macro-Cells as input and output.
 * All configuration is made using an XML-config file and does not require recompilation when modified.
 * @Author Felix Maurer
 */
namespace coupling{
    template<unsigned int dim>
    class FilterPipeline;

	//TODO: comment
	enum class Scope { perInstance, postMultiInstance};
}


template<unsigned int dim>
class coupling::FilterPipeline{
    public:

		//TODO: comment! difference: whole domain vs only md2macro incl. indexing
        FilterPipeline(
			std::vector<coupling::datastructures::MacroscopicCell<dim>* > inputCells,
			const tarch::utils::MultiMDService<dim>& multiMDService,
			const coupling::Scope scope,
			const char* cfgpath);
               
        ~FilterPipeline() {
            for(auto sequence : _sequences) delete sequence;
			//TODO: do i have to delete the _...cells as well?

            #ifdef DEBUG_FILTER_PIPELINE
            std::cout << "FP: FilterPipeline deconstructed." << std::endl;
            #endif
        }


		/*
		 * Applies each FilterSequence in order of their appearance in the config file.
		 * Ouput of the specified output-FilterSequence will be written to _md2MacroCells.
		 */
        void operator()();
		

		/*
		 * Getters for FilterSequences.
		 * Not that Junction is a subtype of Sequence, so this is how to get Junctions as well.
		 */
       	coupling::FilterSequence<dim> * getSequence(const char* identifier) const;
       	std::vector<coupling::FilterSequence<dim> *> getAllSequences() const { return _sequences; }

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
       	void loadSequencesFromXML(tinyxml2::XMLElement* metaNode);

		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _md2MacroCells;
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outerCells;

		const tarch::utils::MultiMDService<dim>& _multiMDService;

		tinyxml2::XMLDocument _config;

		const coupling::Scope _scope;

       	std::vector<coupling::FilterSequence<dim> *> _sequences; 
};


//include implementation of header
#include "FilterPipeline.cpph"
