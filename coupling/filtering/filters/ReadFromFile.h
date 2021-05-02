// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"

#include <string>
#include <vector>
#include <fstream>

#define DEBUG_READ_FROM_FILE

namespace coupling {
    template<unsigned int dim>
    class ReadFromFile;
}


/*
 * Filter that reads cell data from a specified file in .csv format and then writes that data to its output cells.
 * @author Felix Maurer
 *
 * Input format must be compliant to the usual MaMiCo CSV format (using ';' as separator).
 * The following order is assumed:
 * - current iteration
 * - md-to-macro indexing
 * - mamico indexing
 * - scalars
 * - vectors
 *
 * The input file must contain one separate line per cell.
 */
template<unsigned int dim>
class coupling::ReadFromFile : public coupling::FilterInterface<dim>{
    public:
        ReadFromFile(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
				const std::array<bool, 7> filteredValues, 
				std::string location):

				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "RFF"),
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
