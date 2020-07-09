// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>

#define DEBUG_GAUSS
#include "coupling/filtering/FilterInterface.h"

namespace coupling {
    template<unsigned int dim>
    class Gauss;
}


template<unsigned int dim>
class coupling::Gauss : public coupling::FilterInterface<dim>{
    public:
        Gauss(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //covers the entire MD domain
				tarch::la::Vector<dim, unsigned int> domainStart,
				tarch::la::Vector<dim, unsigned int> domainEnd):

				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, domainStart, domainEnd)
		{
        #ifdef DEBUG_GAUSS
			std::cout << "GAUSS: Created Gaussian filter";
        #endif
        }

        ~Gauss(){
        #ifdef DEBUG_GAUSS
            std::cout << "GAUSS: Gaussian filter deconstructed" << std::endl;
        #endif
        }

     
	    void operator()();
};

//include implementation of header
#include "Gauss.cpph"
