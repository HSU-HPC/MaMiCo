// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include "coupling/IndexConversion.h"

//INCLUDE ALL FILTER HEADERS HERE
#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/filters/Gauss.h"
#include "coupling/filtering/filters/POD.h"
#include "coupling/filtering/filters/Strouhal.h"

/*
 * Filter Sequences are used to group filters that will be applied in chronological order.
 * It is possible to customize
 *  - sequence input
 * 	- sequence domain
 * 	- filter parameters
 * per Filter Sequence.
 * @Author Felix Maurer
 */

namespace coupling{
    template<unsigned int dim>
    class FilterSequence;
}

template<unsigned int dim>
class coupling::FilterSequence {
	public:
		/*
		 * Filter Sequences are constructed in coupling::FilterPipeline::loadSequencesFromXML(...).
		 * inputCellVector and cellIndices cover the entire md2Macro domain.
		 * domainStart and domainEnd span a subspace of the global (md2Macro) domain.
		 */
    	FilterSequence( const coupling::IndexConversion<dim>* indexConversion,
						const tarch::utils::MultiMDService<dim>& multiMDService,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	inputCellVector,
						std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
						tarch::la::Vector<dim, unsigned int> domainStart,
						tarch::la::Vector<dim, unsigned int> domainEnd,
						bool filteredValues[7]):
    	_indexConversion(indexConversion), 
		_multiMDService(multiMDService),
		_name(name), 
		_inputCellVector(inputCellVector),
		_cellIndices(cellIndices), 
		_domainStart(domainStart),
		_domainEnd(domainEnd),
		_filteredValues(filteredValues),
		_isOutput(false) //potentially updated via loadSequencesFromXML calling setAsOutput()
		{	
			#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << PRINT_PREFIX() << "Now initializing." << std::endl;
        	#endif
            for(unsigned int d = 0; d < dim; d++) if(_cellIndices.back()[d] < _domainEnd[d]){
                std::cout << "Filter domain size larger than MD domain. Aborting." << std::endl;
                exit(EXIT_FAILURE);
            }
			initCellVectors();
			initDomain();
			#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << PRINT_PREFIX() << "Finished initialization." << std::endl;
        	#endif
    	}


    	~FilterSequence(){
			for (auto v1 : _cellVector1) delete v1;
			for (auto v2 : _cellVector2) delete v2;
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << PRINT_PREFIX() << "Deconstructed." << std::endl;
        	#endif
    	}

		/*
		 * Interprets this sequence's filters' parameters (specified in filter_pipeline.xml) and creates filter objects based on that.
		 * Since different filters require vastly different parameters and thus different means of instanciation,
		 * you need to add each new filter manually to this method's "if/else-chain".
		 *
		 * Since after all filter objects are created it is possible to determine whether _cellVector1 or _cellVector2 will be used as output,
		 * this is also done in here.
		 */
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);
		
		/*
		 * Each sequence operates on their own two copies of the global domain. 
		 * Thus, before applying the sequence, we need to update these two copies.
		 */
		void updateCellVectors(){
			for(unsigned int index = 0; index < _cellIndices.size(); index++){
				*(_cellVector1[index]) = *(_inputCellVector[index]);
				*(_cellVector2[index]) = *(_inputCellVector[index]);
			}
		}

    	const char* getName() { return _name; }

		bool isOutput() { return _isOutput; }
		void setAsOutput() { _isOutput = true; }

		/*
		 * Which one of the two cell vectors are this sequence's output is solely dependant on the number of filters the sequence contains,
		 * because each swap (see loadFiltersFromXML) changes what is the sequence's final filter's output.
		 */
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getOutputCellVector() const{ 
			if(_filters.empty()) std::cout << PRINT_PREFIX() << "Warning: Accessing cell vectors while _filters is empty." << std::endl;
			if(_filters.size() % 2 == 0) return _cellVector1;
			else return _cellVector2;
		}

		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	
      
	private:
		//TODO: move to constructor (?)
		/*
		 * Determines based on _domainStart and _domainEnd which of the global domain's cell belong to the sequence's local domain.
		 * This initializes all domain vector member variables (see below).
		 */
		void initDomain();

		//TODO: move to constructor (?)
		/*
		 * Copies all (global) input cells to _cellVector1 and _cellVector2.
		 */
		void initCellVectors();

		std::string PRINT_PREFIX() const{
			return std::string("	FS(").std::string::append(_name).std::string::append("): ");
		}

		const coupling::IndexConversion<dim>* _indexConversion;
		const tarch::utils::MultiMDService<dim> _multiMDService;

    	const char* _name;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector;//points to (foreign) input vector TODO: remove once init...() functions are part of constructor
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1;//allocated for this sequence only
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector2;//allocated for this sequence only
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;//all of the above use the same indexing

        /*pseudo const*/std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputDomainCellVector;//points to macro cells of this sequence's input that are within domain
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector1;//points to cells of vector 1 that are within the domain's range
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector2;//same for vector 2
		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;
    	std::vector<tarch::la::Vector<dim, unsigned int>> _globalDomainCellIndices; //uses _cellIndices indexing (i. e. _domainStart .... _domainEnd))
    	std::vector<tarch::la::Vector<dim, unsigned int>> _localDomainCellIndices; //starts at (0,...,0)
		
		bool *_filteredValues; //actually, bool[7]

		bool _isOutput; //true if this sequence's output vector (see above) is the Filter Pipeline's output
		
		std::vector<coupling::FilterInterface<dim> *> _filters;
		


};

//inlcude implementation
#include "coupling/filtering/FilterSequence.cpph"
