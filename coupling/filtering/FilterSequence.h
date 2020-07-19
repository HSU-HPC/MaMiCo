// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include "coupling/IndexConversion.h"
#include "tarch/configuration/ParseConfiguration.h"

//INCLUDE ALL FILTERS HERE
#include "coupling/filtering/filters/WriteToFile.h"
#include "coupling/filtering/filters/Gauss.h"

//Filter Sequences are used to logically group filters that will be used in chronological order.
//@Author Felix Maurer

namespace coupling{
    template<unsigned int dim>
    class FilterSequence;
}

template<unsigned int dim>
class coupling::FilterSequence {
	public:
    	FilterSequence( const coupling::IndexConversion<dim>* indexConversion,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	inputCellVector,
						std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
						tarch::la::Vector<dim, unsigned int> domainStart,
						tarch::la::Vector<dim, unsigned int> domainEnd):
    	_indexConversion(indexConversion), 
		_name(name), 
		_inputCellVector(inputCellVector),
		_cellIndices(cellIndices), 
		_domainStart(domainStart),
		_domainEnd(domainEnd),
		_isOutput(false)
		{	
			#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FS (" << _name << "): Now initializing." << std::endl;
        	#endif
            for(unsigned int d = 0; d < dim; d++) if(_cellIndices.back()[d] < _domainEnd[d]){
                std::cout << "Filter domain size larger than MD domain. Aborting." << std::endl;
                exit(EXIT_FAILURE);
            }
			initCellVectors();
			initDomain();
			#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FS ("<< _name <<"): Finished initialization " << std::endl;
        	#endif
    	}


    	~FilterSequence(){
			for (auto v1 : _cellVector1) delete v1;
			for (auto v2 : _cellVector2) delete v2;
			for (auto f : _filters) delete f;
        	#ifdef DEBUG_FILTER_PIPELINE
        	std::cout << "FS("<< _name <<"): Deconstructed. " << std::endl;
        	#endif
    	}

		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode);
					
		void updateCellVectors(){
			for(unsigned int index = 0; index < _cellIndices.size(); index++){
				*(_cellVector1[index]) = *(_inputCellVector[index]);
				*(_cellVector2[index]) = *(_inputCellVector[index]);
			}
		}

    	const char* getName() { return _name; }

		bool isOutput() { return _isOutput; };
		void setAsOutput() { _isOutput = true; };

    	const std::vector<coupling::datastructures::MacroscopicCell<dim>* >& getOutputCellVector() const{ 
			if(_filters.empty()) std::cout << "FS (" << _name << "): Warning: Accessing cell vectors while _filters is empty." << std::endl;
			if(_filters.size() % 2 == 0) return _cellVector1;
			else return _cellVector2;
		}

    	//void setCellIndices(const std::vector<tarch::la::Vector<dim, unsigned int>>& cellIndices) { _cellIndices = cellIndices; }
    	//std::vector<tarch::la::Vector<dim, unsigned int>> getCellIndices() { return _cellIndices; }
       
		std::vector<coupling::FilterInterface<dim> *> getFilters() { return _filters; }	

 		//TODO: remove
		void DEBUG_GET_VECTOR_SIZES(){
			std::cout << "FS (" << _name << "): Printing vector sizes, capacities" << std::endl;
			std::cout << "_cellVector1: " << _cellVector1.size() << " , " << _cellVector1.capacity() << std::endl;
			std::cout << "_cellVector2: " << _cellVector2.size() << " , " << _cellVector2.capacity() << std::endl;
			std::cout << "_cellIndices: " << _cellIndices.size() << " , " << _cellIndices.capacity() << std::endl;
			std::cout << "_inputDomainCellVector: " << _inputDomainCellVector.size() << " , " << _inputDomainCellVector.capacity() << std::endl;
			std::cout << "_domainCellVector1: " << _domainCellVector1.size() << " , " << _domainCellVector1.capacity() << std::endl;
			std::cout << "_domainCellVector2: " << _domainCellVector2.size() << " , " << _domainCellVector1.capacity() << std::endl;
		}       
	private:
		//TODO: move to constructor (?)
		void initDomain();
		//TODO: move to constructor (?)
		void initCellVectors();

		const coupling::IndexConversion<dim>* _indexConversion;
    	const char* _name;

    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector;//points to (foreign) input vector
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1;//allocated for this sequence only
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector2;//allocated for this sequence only
    	std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices;//all of the above use the same indexing

        /*pseudo const*/std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputDomainCellVector;//pointers to macro cells of this sequence's input that are within domain
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector1;//points to cells of vector 1 that are within the domain's range
    	std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCellVector2;//same for vector 2
		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;
    	std::vector<tarch::la::Vector<dim, unsigned int>> _globalDomainCellIndices; //uses _cellIndices indexing (i. e. _domainStart .... _domainEnd))
    	std::vector<tarch::la::Vector<dim, unsigned int>> _localDomainCellIndices; //starts at (0,0,0)
		
		bool _isOutput;
		std::vector<coupling::FilterInterface<dim> *> _filters;
		


};

//inlcude implementation
#include "coupling/filtering/FilterSequence.cpph"
