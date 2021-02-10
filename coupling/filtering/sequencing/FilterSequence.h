// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include<map>
#include<numeric>
#include<sys/time.h>

#include "coupling/IndexConversionMD2Macro.h"
#include "tarch/configuration/ParseConfiguration.h"

//INCLUDE ALL FILTER HEADERS HERE
#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/filters/Gauss.h"
#include "coupling/filtering/filters/POD.h"
#include "coupling/filtering/filters/NLM.h"
#include "coupling/filtering/filters/Strouhal.h"
#include "coupling/filtering/filters/FilterFromFunction.h"
#include "coupling/filtering/filters/Copy.h"
#include "coupling/filtering/SequentialFilter.h"

/*
 * Filter Sequences are used to group filters that will be applied in chronological order.
 * It is possible to customize
 *  - sequence input
 * 	- sequence domain
 * 	- filter parameters
 * 	- modifiablitiy of filter list
 * per Filter Sequence.
 *
 * A generalized version of this concept is FilterJunction.
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
    	FilterSequence( const coupling::IndexConversionMD2Macro<dim>* indexConversion,
						const tarch::utils::MultiMDService<dim>& multiMDService,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	inputCellVector,
						std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
						tarch::la::Vector<dim, unsigned int> domainStart,
						tarch::la::Vector<dim, unsigned int> domainEnd,
						std::array<bool, 7> filteredValues):
    	_ic(indexConversion), 
		_multiMDService(multiMDService),
		_name(name), 
		_inputCellVector(inputCellVector),
		_cellIndices(cellIndices), 
		_domainStart(domainStart),
		_domainEnd(domainEnd),
		_filteredValues(filteredValues),
		_isOutput(false), //potentially updated via loadSequencesFromXML calling setAsOutput()
		_isModifiable(true), //TODO: allow const sequences via XML attribute
		_timestepsElapsed(0)
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

			bool filtersAnything = false;
			for(unsigned int i = 0; i < 7; i++) if(_filteredValues[i]) filtersAnything = true;
			if(!filtersAnything) std::cout << "Warning: Filter sequence " << _name << " does not filter any values. Add 'filtered-values' attribute to XML element to change this." << std::endl;
			#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << PRINT_PREFIX() << "Finished initialization." << std::endl;
        	#endif
    	}


    	virtual ~FilterSequence(){
			for (auto v1 : _cellVector1) delete v1;
			for (auto v2 : _cellVector2) delete v2;
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << PRINT_PREFIX() << "Deconstructed." << std::endl;
        	#endif

			//Output average application times for all filters
			std::cout << PRINT_PREFIX() << "Average application times for filters in this sequence in \u03BCs:" << std::endl;
			for(auto filter : _filters){
				std::cout << "	" << filter->getType() << ": " << (double) std::accumulate(_filterTimes[filter].begin(),_filterTimes[filter].end(),0) / (double) _timestepsElapsed << std::endl;
			}
    	}
	
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

		/*
		 * Applies all filters stored in _filters. Also measures performance times.
		 */
		void operator()() {
			for(auto filter : _filters){
				//TODO: isnt there a C++ way to do this?
				//measure time before application
				timeval before;
				gettimeofday(&before, NULL);

				//Apply the filter's operator()
				(*filter)();

				//measure time after application
				timeval after;
				gettimeofday(&after, NULL);

				//store time difference in usec in map
				_filterTimes[filter].push_back((after.tv_sec*1000000+after.tv_usec) - (before.tv_sec*1000000+before.tv_usec));

				#ifdef DEBUG_FILTER_PIPELINE
					std::cout 
						<< PRINT_PREFIX() << "Applied filter " << filter->getType() 
						<< ". Application time: " << _filterTimes[filter].back() << "\u03BCs" << std::endl;
				#endif
			}
			_timestepsElapsed++;
		}

		/*
		 * GETTER/SETTER SECTION
		 */

    	const char* getName() { return _name; }

		bool isOutput() { return _isOutput; }
		//TODO: is this used?
		void setAsOutput() { 
			std::cout << PRINT_PREFIX() << " Setting as pipeline to macro solver output." << std::endl;
			_isOutput = true; 
		}

		bool isModifiable() { return _isModifiable; }
		void makeUnmodifiable() { _isModifiable = false; }

		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	

		/*
		 * All virtual functions below are redefined in case this sequence is actually a FilterJunction.
		 * Read FilterJunction.h carefully.
		 */

		/*
		 * Interprets this sequence's filters' parameters (specified in filter_pipeline.xml) and creates filter objects based on that.
		 * Since different filters require vastly different parameters and thus different means of instanciation,
		 * you need to add each new filter manually to this method's "if/else-chain".
		 *
		 * Since after all filter objects are created it is possible to determine whether _cellVector1 or _cellVector2 will be used as output,
		 * this is also done in here.
		 *
		 * In addition to that, if this sequence is declared as unmodifibale, this gets also detected in here
		 */
		virtual int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);

		/*
		 * Determines based on _domainStart and _domainEnd which of the global domain's cell belong to the sequence's local domain.
		 * This initializes all domain vector member variables (see below).
		 *
		 * Used in constructor. TODO: i like to move it
		 *
		 * Implemented by FilterJunction to init domain partitions.
		 */
		virtual void initDomain();

		/*
		 * This member function allows appendance and insertion of filters defined by two processing functions to a modifiable sequence at runtime. Index -1 implies appending. 
		 */
		virtual void addFilter( 	
				const std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* applyScalar,
				const std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* applyVector,
				int filterIndex = -1
		);

		/*
		 * Which one of the two cell vectors are this sequence's output is solely dependant on the number of filters the sequence contains,
		 * because each swap (see loadFiltersFromXML) changes what is the sequence's final filter's output.
		 *
		 * Some sequences have more than one output, thus the optional parameter. Has no effect on a basic FilterSequence.
		 */
    	virtual const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getOutputCellVector(unsigned int outputIndex = 0) const{ 
			if(_filters.empty()) std::cout << PRINT_PREFIX() << "Warning: Accessing cell vectors while _filters is empty." << std::endl;
			if(_filters.size() % 2 == 0) return _cellVector1;
			else return _cellVector2;
		}

		virtual void printOutputCellVector() const {
			(getOutputCellVector() == _cellVector1) ? std::cout << "Cell vector 1 " : std::cout << "Cell vector 2 ";
			std::cout << "will be used as output." << std::endl;
		}
		
		/*
		 * This is trivial for the traditional sequence, i.e. non-Junction-case.
		 */
		virtual unsigned int getNumInputs() { return 1; }
		virtual unsigned int getNumOutputs() { return 1; }



		/*
		 * DUMMY HELPER FUNCTIONS
		 */

		virtual void printFilters() {
			std::cout << "Filters in sequence " << _name << ": ";
			for(auto f : _filters) std::cout << f->getType() << " ";
			std::cout << std::endl;
		}

		void pintOutputVector() const {
			std::cout << PRINT_PREFIX() << "Number of Filters: " << _filters.size() <<". Output vector will be ";
			if(_filters.size() % 2 == 0)  std::cout << "_cellVector1." << std::endl;
			else std::cout << "_cellVector2." << std::endl;
		}

		virtual std::string PRINT_PREFIX() const {
			return std::string("	FS(").std::string::append(_name).std::string::append("): ");
		}

      
	private:
		/*
		 * Copies all (global) input cells to _cellVector1 and _cellVector2.
		 *
		 * Used in consctructor. //TODO: move it
		 *
		 * Not implemented by FilterJunction: Equivalent procedure can be found in that class' constructor.
		 */
		void initCellVectors();
	
	protected:

		const coupling::IndexConversionMD2Macro<dim>* _ic;
		const tarch::utils::MultiMDService<dim> _multiMDService;

    	const char* _name;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector;//points to (foreign) input vector 
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1;//allocated for this sequence only
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector2;//allocated for this sequence only
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;//all of the above use the same indexing

		//points to macro cells of this sequence's input that are within domain
        /*pseudo const*/std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputDomainCellVector;
		//points to cells of vector 1 that are within the domain's range
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector1;
		//same for vector 2
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector2;

		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;

		//uses _cellIndices indexing (i. e. _domainStart .... _domainEnd))
    	std::vector<tarch::la::Vector<dim, unsigned int>> _globalDomainCellIndices;     	
		//starts at (0,...,0)
		std::vector<tarch::la::Vector<dim, unsigned int>> _localDomainCellIndices; 
		
		std::array<bool, 7> _filteredValues;

		bool _isOutput; //true if this sequence's output vector (see above) is the Filter Pipeline's output
		bool _isModifiable; //true while filters can be added to sequence
		
		std::vector<coupling::FilterInterface<dim> *> _filters;

		//stores application times for all filters
		std::map<coupling::FilterInterface<dim> *, std::vector<unsigned int>> _filterTimes;

		//there must be some other place to get this from
		unsigned int _timestepsElapsed;
};

//inlcude implementation
#include "coupling/filtering/sequencing/FilterSequence.cpph"
