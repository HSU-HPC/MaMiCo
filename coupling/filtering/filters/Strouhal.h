// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/FilterInterfaceReadOnly.h"

#include <vector>

//#define DEBUG_STROUHAL

namespace coupling {
    template<unsigned int dim>
    class Strouhal;
}

/**
 * TODO: Comment
 */
template<unsigned int dim>
class coupling::Strouhal : public coupling::FilterInterfaceReadOnly<dim>{
    public:
        Strouhal(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
				bool filteredValues[7],
				double u,
				double d):
				coupling::FilterInterfaceReadOnly<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues),
				_U(u),
				_D(d)
		{	
			if(dim < 2){
				std::cout << "ERROR: Strouhal filter only works for dim >= 2." << std::endl;
				exit(EXIT_FAILURE);
			}
        	#ifdef DEBUG_STROUHAL
            std::cout << "		STROUHAL: Instance created." << std::endl;
        	#endif
        }

        ~Strouhal(){
			std::cout << "STROUHAL NUMBER IN MD DOMAIN: " << calculateStrouhalNumber() << std::endl;
        	#ifdef DEBUG_WRITE_TO_FILE
            std::cout << "		STROUHAL: Instance deconstructed." << std::endl;
        	#endif
        }

     
	    void operator()();

    private:
		double calculateStrouhalNumber();

		std::vector<double> _v_y;
		double _U;
		double _D;
};

//include implementation of header
#include "Strouhal.cpph"
