// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_JUNCTION

#include "coupling/filtering/sequencing/FilterSequence.h"

//INCLUDE ALL JUNCTOR HEADERS HERE
#include "coupling/filtering/filters/NLM.h"
#include "coupling/filtering/filters/AveragingJunctor.h"
#include "coupling/filtering/interfaces/JunctorInterface.h" //this will be redundant in the future

/*
 * WORK IN PROGRESS. USE WITH CAUTION
 *
 * Generalizes the concept of FilterSequences: A FilterJunction can have more than a single input.
 * This allows for seemless integration of FilterJunctions into a net of FilterSequences.
 *
 * TODO: 
 *  - Support multiple outputs. ("X-Junctions")
 *  - Support dynamically linked filters.
 * @Author Felix Maurer
 */

namespace coupling{
	template<unsigned int dim, std::size_t inputc>
    class FilterJunction;
}

template<unsigned int dim, std::size_t inputc>
class coupling::FilterJunction : public coupling::FilterSequence<dim> {
	public:
    	FilterJunction( const coupling::IndexConversionMD2Macro<dim>* indexConversion,
						const tarch::utils::MultiMDService<dim>& multiMDService,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	inputCellVector, //concatenation of numberImput input cell vectors
						std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, 
						tarch::la::Vector<dim, unsigned int> domainStart,
						tarch::la::Vector<dim, unsigned int> domainEnd,
						std::array<bool, 7> filteredValues
		):
		coupling::FilterSequence<dim>(indexConversion, multiMDService, name, inputCellVector, cellIndices, domainStart, domainEnd, filteredValues)
		{	
			if(inputc == 0)
				throw std::runtime_error("ERROR: Creating FilterJunction with inputc = 0.");

			#ifdef DEBUG_FILTER_JUNCTION
        	std::cout << PRINT_PREFIX() << "This is a FilterJunction. Number of inputs:" << inputc << std::endl;
        	#endif

			//Partition input vector
			unsigned int partitionSize = coupling::FilterSequence<dim>::_inputCellVector.size() / inputc;

			for( unsigned int p = 0; p < inputc; p++) {
				//_inputCellVector 
				_inputCellVector_parted[p] =
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_inputCellVector.begin() + (p * partitionSize), 
					coupling::FilterSequence<dim>::_inputCellVector.begin() + ((p+1) * partitionSize));
				#ifdef DEBUG_FILTER_JUNCTION
					std::cout << PRINT_PREFIX() << "Size of _inputCellVector_parted[" << p << "]: " << _inputCellVector_parted[p].size() << std::endl;
				#endif


				//_cellVector1
				_cellVector1_parted[p] =
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_cellVector1.begin() + (p * partitionSize), 
					coupling::FilterSequence<dim>::_cellVector1.begin() + ((p+1) * partitionSize));
				#ifdef DEBUG_FILTER_JUNCTION
					std::cout << PRINT_PREFIX() << "Size of _cellVector1_parted[" << p << "]: " << _cellVector1_parted[p].size() << std::endl;
				#endif

				//_cellVector2
				_cellVector2_parted[p] =
					std::vector<coupling::datastructures::MacroscopicCell<dim> *>
					(coupling::FilterSequence<dim>::_cellVector2.begin() + (p * partitionSize), 
					coupling::FilterSequence<dim>::_cellVector2.begin() + ((p+1) * partitionSize));
				#ifdef DEBUG_FILTER_JUNCTION
					std::cout << PRINT_PREFIX() << "Size of _cellVector2_parted[" << p << "]: " << _cellVector2_parted[p].size() << std::endl;
				#endif
			}

			initDomain();

			coupling::FilterSequence<dim>::_isModifiable = false; //Dynamic filters are not yet supported.
    	}

    	~FilterJunction(){}

		/*
		 * This member function allows appendance and insertion of filters defined by two processing functions to a modifiable sequence at runtime. Index -1 implies appending. 
		 */
		void addFilter( 	
				const std::function<std::vector<double> (std::vector<double>, std::vector<std::array<unsigned int, dim>>)>* applyScalar,
				const std::function<std::vector<std::array<double, dim>> (std::vector<std::array<double, dim>>, std::vector<std::array<unsigned int, dim>>)>* applyVector,
				int filterIndex = -1
		) override {
			//Do nothing, not yet supported.
			#ifdef DEBUG_FILTER_JUNCTION
			std::cout << PRINT_PREFIX() << "This is a FilterJunction. addFilter(...) is not supported and has no effect." << std::endl;
			#endif
		}

		/*
		 * This function is very similar to the interface's. Check coupling::FilterSequence for more details.
		 */
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode) override;

		/*
		 * This function is very similar to the interface's. Check coupling::FilterSequence for more details.
		 */
		void initDomain() override;

		/*
		 * The first partition of _cellVector1/2 is the main partition. A junction's default output is always its main partition.
		 */
    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getOutputCellVector(unsigned int outputIndex = 0) const override{ 
			if(outputIndex >= inputc) {
				std::cout << PRINT_PREFIX() << "ERROR: getOutputCellVector: Requested output index(" << outputIndex << ") too high. (partitions: )" << inputc << std::endl;
				exit(EXIT_FAILURE);
			}

			if(coupling::FilterSequence<dim>::_filters.empty()) std::cout << PRINT_PREFIX() << "Warning: Accessing cell vectors while _filters is empty." << std::endl;
			if(coupling::FilterSequence<dim>::_filters.size() % 2 == 0) return _cellVector1_parted[0];
			else return _cellVector2_parted[0];
		}


		void printFilters() override {
			std::cout << "Junctors in junction " << coupling::FilterSequence<dim>::_name << ": ";
			for(auto f : coupling::FilterSequence<dim>::_filters) std::cout << f->getType() << " ";
			std::cout << std::endl;
		}

		std::string PRINT_PREFIX() const override {
			return std::string("	FJ(").std::string::append(coupling::FilterSequence<dim>::_name).std::string::append("): ");
		}


	private:
		//These must be parted for junction output and junctor in/output.
		std::array<std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _inputCellVector_parted;
		std::array<std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _cellVector1_parted;
		std::array<std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _cellVector2_parted;

		//Partitions are given to junctors. Goal: Arbitrary number of input/output partitions.
		std::array</*pseudo const*/std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _inputDomainCellVector_parted;
		std::array<std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _domainCellVector1_parted;
		std::array<std::vector<coupling::datastructures::MacroscopicCell<dim>* >, inputc> _domainCellVector2_parted;
};

//inlcude implementation
#include "coupling/filtering/sequencing/FilterJunction.cpph"
