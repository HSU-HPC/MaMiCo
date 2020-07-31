// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>

//#define DEBUG_GAUSS
#include "coupling/filtering/FilterInterface.h"

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

template<unsigned int dim>
class coupling::Gauss : public coupling::FilterInterface<dim>{
    public:
        Gauss(  const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //Use local indexing! (starting at (0,...,0))
				unsigned int d,
				const bool uses[7]):
				
				
				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices),
				_dim(d),
				_lastIndex(coupling::FilterInterface<dim>::_cellIndices.back()),
				_useMicroMass(uses[0]),
			 	_useMicroMomentum(uses[1]),
				_useMacroMass(uses[2]),
				_useMacroMomentum(uses[3]),
				_usePotEnergy(uses[4]),
				_useVelocity(uses[5]),
				_useTemperature(uses[6])
		{
			for(unsigned int d = 0; d < dim; d++)if(coupling::FilterInterface<dim>::_cellIndices.back()[d] < 2){
				std::cout << "ERROR: GAUSS: Input domain has to be " << dim << " dimensional." << std::endl;
			}
        	#ifdef DEBUG_GAUSS
			std::cout << "		GAUSS: Created Gaussian filter. It will filter:" << std::endl;
			if(_useMicroMass) std::cout << "		- Microscopic mass" << std::endl;
			if(_useMicroMomentum) std::cout << "		- Microscopic momentum" << std::endl;
			if(_useMacroMass) std::cout << "		- Macroscopic mass" << std::endl;
			if(_useMacroMomentum) std::cout << "		- Macroscopic momentum" << std::endl;
			if(_usePotEnergy) std::cout << "		- Potential energy" << std::endl;
			if(_useVelocity) std::cout << "		- Velocity" << std::endl;
			if(_useTemperature) std::cout << "		- Temperature" << std::endl;
       		#endif
        }

        ~Gauss(){
        	#ifdef DEBUG_GAUSS
            std::cout << "		GAUSS: Gaussian filter deconstructed" << std::endl;
        	#endif
        }

     
	    void operator()();
	private:
		//on which axis this filter operates. 0 <= _dim <= dim
		unsigned int _dim;

		//returns the cell that's above the cell at index on the d-axis
		unsigned int getLowerCellIndex(unsigned int index, unsigned int d);
		//returns the cell that's below the cell at index on the d-axis
		unsigned int getUpperCellIndex(unsigned int index, unsigned int d);

		//TODO: Move to interface?
		tarch::la::Vector<dim, unsigned int> _lastIndex;

		//TODO: Move to interface?
		bool _useMicroMass;
		bool _useMicroMomentum;
		bool _useMacroMass;
		bool _useMacroMomentum;
		bool _usePotEnergy;
		bool _useVelocity;
		bool _useTemperature;
};

//include implementation of header
#include "Gauss.cpph"
