// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>
#include <cmath>

//#define DEBUG_GAUSS
#include "coupling/filtering/interfaces/FilterInterface.h"

namespace coupling {
    template<unsigned int dim>
    class Gauss;

	//cf. member variable in coupling::Gauss for more details
	enum GaussExtrapolationStrategy {NONE, LINEAR};
}

//Define kernel radius. e.g. radius = 1 means kernel size of 3
#define GAUSS_KERNEL_RADIUS 1

/*
 * Implements a gaussian filter. Limited to what is listed in the disclaimer above. 
 * Only specified in one dimension: If you wish to use a multidimensional gaussian filter, simply chain multiple instances of this filter in one FilterSequence.
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::Gauss : public coupling::FilterInterface<dim>{
	using coupling::FilterInterface<dim>::_inputCells;
	using coupling::FilterInterface<dim>::_outputCells;
	using coupling::FilterInterface<dim>::_cellIndices;
	using coupling::FilterInterface<dim>::_scalarGetters;
	using coupling::FilterInterface<dim>::_vectorGetters;
	using coupling::FilterInterface<dim>::_scalarSetters;
	using coupling::FilterInterface<dim>::_vectorSetters;

    public:
        Gauss(  const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //Use local indexing! (starting at (0,...,0))
				const std::array<bool, 7> filteredValues,
				unsigned int dimension,
				unsigned int sigma,
				const char* extrapolationStrategy):
				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "GAUSS"),
				_dim(dimension),
				_sigma(sigma),
				_lastIndex(coupling::FilterInterface<dim>::_cellIndices.back()),
				_kernel(generateKernel())
		{
			std::cout << "WARNING: You're using a GAUSS-Filter. As this filter has not been tested thoroughly, caution is advised!" << std::endl;

			//TODO
			if(GAUSS_KERNEL_RADIUS != 1) 
				throw std::runtime_error("ERROR: GAUSS: Kernel radius != 1 currently not supported.");

			//TODO
			if(sigma != 1) 
				throw std::runtime_error("ERROR: GAUSS: sigma != 1 currently not supported.");

			//Overwrite kernel with hardcoded weights TODO: fix kernel generation, then remove this
			_kernel = {0.27901, 0.44198, 0.27901};


			if(coupling::FilterInterface<dim>::_cellIndices.back()[_dim] < 2)
				throw std::runtime_error("ERROR: GAUSS: Invalid input domain.");

			if(extrapolationStrategy == nullptr || std::strcmp(extrapolationStrategy, "none") == 0) _extrapolationStrategy = NONE;
			else if(std::strcmp(extrapolationStrategy, "linear") == 0) _extrapolationStrategy = LINEAR;
			else {
				std::cout << "Extrapolation strategy: " << extrapolationStrategy << std::endl;
				throw std::runtime_error("ERROR: GAUSS: Unknown extrapolation strategy.");
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
		std::array<double, 1+2*GAUSS_KERNEL_RADIUS> generateKernel();

		double gaussianDensityFunction(int x);

		//returns the cell that's above the cell at index on the d-axis
		unsigned int getIndexBelow(unsigned int index, unsigned int d);
		//returns the cell that's below the cell at index on the d-axis
		unsigned int getIndexAbove(unsigned int index, unsigned int d);


		//on which axis this filter operates. 0 <= _dim <= dim
		unsigned int _dim;

		//standard deviation used
		unsigned int _sigma;

		tarch::la::Vector<dim, unsigned int> _lastIndex;

		std::array<double, 1+2*GAUSS_KERNEL_RADIUS> _kernel;
		
		/**
		 * Determines how to apply filter to border cells:
		 * NONE = only use existing cells and increase their weight accordingly
		 * LINEAR = linear extrapolation
		 */
		coupling::GaussExtrapolationStrategy _extrapolationStrategy;


};

//include implementation of header
#include "Gauss.cpph"
