// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"

#include <string>
#include <vector>
#include <fstream>

//#define DEBUG_READ_FROM_FILE

namespace coupling {
    namespace filtering {
        template<unsigned int dim>
        class ReadFromFile;
    }
}


/*
 * Filter that reads cell data from a specified file in .csv format and then writes that data to its output cells.
 *
 * Input format must be compliant to the usual MaMiCo CSV format (using ';' as separator).
 * The following order is assumed:
 * - current iteration
 * - scalar cell properties
 * - vector cell properties
 *
 * The input file must contain one separate line per cell.
 *
 * @Author Felix Maurer
 */
template<unsigned int dim>
class coupling::filtering::ReadFromFile : public coupling::filtering::FilterInterface<dim>{
    public:
        ReadFromFile(
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
			const std::array<bool, 7> filteredValues, 
			std::string location
		):
			coupling::filtering::FilterInterface<dim>(inputCellVector, outputCellVector, filteredValues, "RFF"),
		   	_location(location),
			_iteration(0)
		{	
        	#ifdef DEBUG_READ_FROM_FILE
            std::cout << "		RFF: Read from file instance created. Will read from: " << _location << "." << std::endl;
        	#endif
        }

        ~ReadFromFile(){
        	#ifdef DEBUG_READ_FROM_FILE
            std::cout << "		RFF: Read from file instance deconstructed." << std::endl;
        	#endif
        }
     
	    void operator()();

    private:
        std::string _location;
		unsigned int _iteration;

        std::ifstream _file;
};

//include implementation of header
#include "ReadFromFile.cpph"
