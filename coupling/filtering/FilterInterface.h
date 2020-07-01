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
		/*FilterInterface(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> inputCellIndices): //covers the entire MD domain
				_domainCells(inputCellVector),
				_domainCellIndices(inputCellIndices)//,
				_domainStart(inputCellIndices[0]),
				_domainEnd(inputCellIndices.back()),
				_domainSize(_domainEnd - _domainStart)
		{}*/

		//Use this in case your filter uses a subdomain of the input domain
		FilterInterface(
				tarch::la::Vector<dim, unsigned int> domainStart, //unit: number of macroscopic cells
				tarch::la::Vector<dim, unsigned int> domainEnd):
				_domainStart(domainStart),
				_domainEnd(domainEnd),
				_domainSize(_domainEnd - _domainStart),
				_isInitialized(false)
		{}

		virtual ~FilterInterface(){};

		//Initializes domain vectors. "inputCellVector/-Indices" must cover the entire MD domain.
		void initDomain(
					const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
					const std::vector<tarch::la::Vector<dim, unsigned int>> inputCellIndices
		){
			#ifdef DEBUG_FILTER
			std::cout << "Filter: Initizing filter. Input vector has " << inputCellVector.size() << " cells. ";
			#endif
			if(inputCellVector.size() != inputCellIndices.size()){
				std::cout << "Cell and index vector out of synch. Aborting.";
				exit(EXIT_FAILURE);
			}
			
			for(unsigned int d = 0; d < dim; d++) if(inputCellIndices.back()[d] < _domainEnd[d]){
				std::cout << "Filter domain size larger than MD domain. Aborting." << std::endl;
				exit(EXIT_FAILURE);
			}

			bool outOfBounds;
			for(unsigned int index = 0; index < inputCellVector.size(); index++){
				outOfBounds = false;
				for(unsigned int d = 0; d < dim; d++) 
					//"index" does not reference the dim-dimensional indices stored in e.g. inputCellIndices but the inputCells' linear index.
					if(inputCellIndices[index][d] < _domainStart[d] || inputCellIndices[index][d] > _domainEnd[d]) outOfBounds = true;
				if(!outOfBounds){
					_domainCells.push_back(inputCellVector[index]);
					_domainCellIndices.push_back(inputCellIndices[index]);
				}
			}//index			

			_isInitialized = true;

			#ifdef DEBUG_FILTER
			std::cout << "Domain vector has " << _domainCells.size() << " cells. " << std::endl;
			#endif

		}

		//Applies the filter to all cells that are within the filter's domain.
		virtual void apply() = 0;

		//TODO: As soon as loadSequenceData() is only called once (and not during FilterPipeline::apply()), we can move this boolean to FilterSequence
		//Checks if the "initDomain(...)" function has been called before and thus can be skipped
		bool isInitialized() { return _isInitialized; }

	protected:
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _domainCells; //use these as input cells for filtering
		std::vector<tarch::la::Vector<dim,unsigned int>> _domainCellIndices;


		//filter domain is determined by spanning a <dim>-dimensional space between two vectors.
		//initiliazed within the filter's constructor, uses XML attributes of the corresponding filter node as input
		tarch::la::Vector<dim, unsigned int> _domainStart;
		tarch::la::Vector<dim, unsigned int> _domainEnd;
		tarch::la::Vector<dim, unsigned int> _domainSize;

	private:
		bool _isInitialized;
};
