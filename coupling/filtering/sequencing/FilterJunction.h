// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_JUNCTION

#include "coupling/filtering/sequencing/FilterSequence.h"

//INCLUDE ALL JUNCTOR HEADERS HERE
#include "coupling/filtering/filters/NLM.h"
#include "coupling/filtering/interfaces/JunctorInterface.h" //this will be redundant in the future

/*
 * WORK IN PROGRESS. USE WITH CAUTION
 *
 * TODO: very big comment
 *
 *
 * TODO: 
 *  - Support multiple outputs. ("X-Junctions")
 *  - Support dynamically linked filters.
 *  - Support custom selection i/o partitions per junctor.
 * @Author Felix Maurer
 */

namespace coupling{
	template<unsigned int dim, std::size_t inputc, std::size_t outputc>
    class FilterJunction;
}

template<unsigned int dim, std::size_t inputc, std::size_t outputc>
class coupling::FilterJunction : public coupling::FilterSequence<dim> {
	public:
    	FilterJunction( const coupling::IndexConversionMD2Macro<dim>* indexConversion,
						const tarch::utils::MultiMDService<dim>& multiMDService,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	inputCellVector, //concatenation of numberImput input cell vectors
						std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //dont concat indices! TODO: check this in constructor
						tarch::la::Vector<dim, unsigned int> domainStart,
						tarch::la::Vector<dim, unsigned int> domainEnd,
						std::array<bool, 7> filteredValues
		):
		coupling::FilterSequence<dim>(indexConversion, multiMDService, name, inputCellVector, cellIndices, domainStart, domainEnd, filteredValues)
		{	
			#ifdef DEBUG_FILTER_JUNCTION
        	std::cout << coupling::FilterSequence<dim>::PRINT_PREFIX() << "This is a FilterJunction. Number of inputs:" << inputc << std::endl;
        	#endif

			//Partition input vector
			//TODO: testing (a lot of it)
			
			//Used for partitioning domain vectors. TODO: check for redundancy with some field of FilterSequence (?)
			unsigned int totalSize = coupling::FilterSequence<dim>::_inputDomainCellVector.size() / inputc;
			unsigned int domainSize = coupling::FilterSequence<dim>::_inputDomainCellVector.size() / inputc;

			for( unsigned int p = 0; p < inputc; p++) {
				//_inputCellVector 
				//TODO Do i need this?

				//_cellVector1
				_cellVector1_parted[p] = (
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_cellVector1.begin() + (p * totalSize)), 
					(coupling::FilterSequence<dim>::_cellVector1.begin() + (p+1 * totalSize))
				);

				//_cellVector2
				_cellVector2_parted.push_back(
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_cellVector2.begin() + (p * totalSize)), 
					(coupling::FilterSequence<dim>::_cellVector2.begin() + (p+1 * totalSize))
				);

				//_inputDomainCellVector
				_inputDomainCellVector_parted.push_back(
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_inputDomainCellVector.begin() + (p * domainSize)), 
					(coupling::FilterSequence<dim>::_inputDomainCellVector.begin() + (p+1 * domainSize))
				);

				//_domaincellVector1
				_domainCellVector1_parted.push_back(
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_domainCellVector1.begin() + (p * domainSize)), 
					(coupling::FilterSequence<dim>::_domainCellVector1.begin() + (p+1 * domainSize))
				);

				//_domainCellVector2
				_domainCellVector2_parted.push_back(
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_domainCellVector2.begin() + (p * domainSize)), 
					(coupling::FilterSequence<dim>::_domainCellVector2.begin() + (p+1 * domainSize))
				);
			}

			coupling::FilterSequence<dim>::_isModifiable = false; //Dynamic filters are not yet supported. TODO
    	}

    	~FilterJunction(){}

		/*
		 * This member function allows appendance and insertion of filters defined by two processing functions to a modifiable sequence at runtime. Index -1 implies appending. 
		 */
		void addFilter( 	
				const std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* applyScalar,
				const std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* applyVector,
				int filterIndex = -1
		){
			//Do nothing, not yet supported. TODO
			#ifdef DEBUG_FILTER_JUNCTION
			std::cout << coupling::FilterSequence<dim>::PRINT_PREFIX() << "This is a FilterJunction. addFilter(...) is not supported and has no effect." << std::endl;
			#endif
		}

		/*
		 * This function is very similar to the interface's. Check coupling::FilterSequence for more details.
		 */
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);

		/*
		 * The first partition of _cellVector1/2 is the main partition. A junction's output is always its main partition. TODO: Multiple outputs.
		 */
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getOutputCellVector() const{ 
			if(coupling::FilterSequence<dim>::_filters.empty()) std::cout << coupling::FilterSequence<dim>::PRINT_PREFIX() << "Warning: Accessing cell vectors while _filters is empty." << std::endl;
			if(coupling::FilterSequence<dim>::_filters.size() % 2 == 0) return _cellVector1_parted[0];
			else return _cellVector2_parted[0];
		}

		std::string PRINT_PREFIX() const {
			return std::string("	FJ(").std::string::append(coupling::FilterInterface<dim>::_name).std::string::append("): ");
		}


	private:

		//TODO: maybe use arrays instead of std::vectors for partitioning?

		//Do I ever need instances of this?
		std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim>* >> _inputCellVector_parted;
		//These must be parted for junction output.
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1_parted[inputc];
		std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim>* >> _cellVector2_parted;

		//Partitions are given to junctors. Goal: Arbitrary number of input/output partitions.
		std::vector</*pseudo const*/std::vector<coupling::datastructures::MacroscopicCell<dim>* >> _inputDomainCellVector_parted;
		std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim>* >> _domainCellVector1_parted;
		std::vector<std::vector<coupling::datastructures::MacroscopicCell<dim>* >> _domainCellVector2_parted;
};

//inlcude implementation
#include "coupling/filtering/sequencing/FilterJunction.cpph"
