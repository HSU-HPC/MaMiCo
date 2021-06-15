// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>

//#define DEBUG_GAUSS
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
    template<unsigned int dim>
    class Gauss;
}

//DISCLAIMER:
//This filter always uses the following parameters for kernel generation:
#define GAUSS_SIGMA 1
#define GAUSS_KERNEL_SIZE 3
//Which implies:
#define GAUSS_OUTER_WEIGHT 0.27901
#define GAUSS_INNER_WEIGHT (1- (GAUSS_OUTER_WEIGHT * 2))

/*
 * Implements a gaussian filter. Limited to what is listed in the disclaimer above. 
 * Only specified in one dimension: If you wish to use a multidimensional gaussian filter, simply chain multiple instances of this filter in one FilterSequence.
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::Gauss : public coupling::FilterInterface<dim>{
    public:
        Gauss(  const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //Use local indexing! (starting at (0,...,0))
				const std::array<bool, 7> filteredValues,
				unsigned int dimension,
				const char* extrapolationStrategy):
				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "GAUSS"),
				_dim(dimension),
				_lastIndex(coupling::FilterInterface<dim>::_cellIndices.back())
		{
			//TODO
			std::cout << "WARNING: You're using a GAUSS-Filter. As this filter has not been tested thoroughly, caution is advised!" << std::endl;

			if(coupling::FilterInterface<dim>::_cellIndices.back()[_dim] < 2){
				std::cout << "ERROR: GAUSS: Invalid input domain: " <<_dim<< std::endl;
				exit(EXIT_FAILURE);
			}

			if(extrapolationStrategy == nullptr || std::strcmp(extrapolationStrategy, "none") == 0) _extrapolationStrategy = 0;
			else if(std::strcmp(extrapolationStrategy, "linear") == 0) _extrapolationStrategy = 1;
			else{
				#ifdef DEBUG_GAUSS
				std::cout << "ERROR: GAUSS: Unknown extrapolation strategy:" << extrapolationStrategy << std::endl;
				#endif
				exit(EXIT_FAILURE);
			}

        	#ifdef DEBUG_GAUSS
			std::cout << "		GAUSS (Dim: " << _dim << "): Created Gaussian filter." << std::endl;
			if(_extrapolationStrategy == 0) std::cout << "		It will not use extrapolation." << std::endl;
			if(_extrapolationStrategy == 1) std::cout << "		It will use linear extrapolation." << std::endl;
       		#endif
        }

        ~Gauss(){
        	#ifdef DEBUG_GAUSS
            std::cout << "		GAUSS (Dim: " << _dim << "): Gaussian filter deconstructed" << std::endl;
        	#endif
        }

     
	    void operator()();
	private:
		//on which axis this filter operates. 0 <= _dim <= dim
		unsigned int _dim;
		
		/**
		 * Determines how to apply filter to border cells:
		 * 0 = only use existing cells and increase their weight accordingly
		 * 1 = linear extrapolation
		 */
		unsigned int _extrapolationStrategy;

		//returns the cell that's above the cell at index on the d-axis
		unsigned int getIndexBelow(unsigned int index, unsigned int d);
		//returns the cell that's below the cell at index on the d-axis
		unsigned int getIndexAbove(unsigned int index, unsigned int d);

		tarch::la::Vector<dim, unsigned int> _lastIndex;
};

//include implementation of header
#include "Gauss.cpph"
