// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER

namespace coupling{
    template<unsigned int dim>
    class FilterInterface;
}


template<unsigned int dim>
class coupling::FilterInterface{
	public:

		//Use this in case your filter uses the entire input domain
		//TODO: Problem: The hosting FilterSequence() does not know its cell vector (inputCellVector) or its size during initialization
		FilterInterface(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices): //covers the entire MD domain
				
				_domainInputCells(inputCellVector),
				_domainOutputCells(outputCellVector),
				_domainCellIndices(cellIndices),
				_domainStart(cellIndices[0]),
				_domainEnd(cellIndices.back())
		{}

		//Use this in case your filter uses a subdomain of the input domain
		FilterInterface(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //covers the entire MD domain
				tarch::la::Vector<dim, unsigned int> domainStart, //unit: number of macroscopic cells
				tarch::la::Vector<dim, unsigned int> domainEnd):
				
				_domainStart(domainStart),
				_domainEnd(domainEnd)
		{
			initDomain(inputCellVector, outputCellVector, cellIndices);
		}

		virtual ~FilterInterface(){};

		
		//Applies the filter to all cells that are within the filter's domain.
		virtual void operator()() = 0;
	protected:
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainInputCells; //pointers to (parts of) either seqCellVec1 or seqCellVec2	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainOutputCells; //pointers to (parts of) of the other seqCellVec
		std::vector<tarch::la::Vector<dim,unsigned int>> _domainCellIndices; //indices for both domainCells vectors is the same


		//filter domain is determined by spanning a <dim>-dimensional space between two vectors.
		//initiliazed within the filter's constructor, uses XML attributes of the corresponding filter node as input
		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;

	private:
		void initDomain(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices
		){
			#ifdef DEBUG_FILTER
			std::cout << "FilterInterface: Initializing filter. Input vector has " << inputCellVector.size() << " cells. ";
			#endif
			//sanity check
			if(inputCellVector.size() != cellIndices.size() || outputCellVector.size() != cellIndices.size()){
				std::cout << "Cell and index vector out of synch. Aborting.";
				exit(EXIT_FAILURE);
			}
		
			//you could just as well use outputCellVector instead of inputCellIndices to iterate	
			for(unsigned int d = 0; d < dim; d++) if(cellIndices.back()[d] < _domainEnd[d]){
				std::cout << "Filter domain size larger than MD domain. Aborting." << std::endl;
				exit(EXIT_FAILURE);
			}
			bool outOfBounds;
			for(unsigned int index = 0; index < inputCellVector.size(); index++){
				outOfBounds = false;
				for(unsigned int d = 0; d < dim; d++) 
					//"index" does not reference the dim-dimensional indices stored in e.g. cellIndices but the cell vectors' linear index.
					if(cellIndices[index][d] < _domainStart[d] || cellIndices[index][d] > _domainEnd[d]) outOfBounds = true;
				if(!outOfBounds){
					_domainInputCells.push_back(inputCellVector[index]);
					_domainOutputCells.push_back(outputCellVector[index]);
					_domainCellIndices.push_back(cellIndices[index]);
				}
			}

			#ifdef DEBUG_FILTER
			std::cout << "Domain vector has " << _domainInputCells.size() << " cells." << std::endl;
			#endif

		}

};
