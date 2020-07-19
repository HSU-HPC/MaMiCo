// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>
#include <fstream>

#define DEBUG_WRITE_TO_FILE
#include "coupling/filtering/FilterInterface.h"

namespace coupling {
    template<unsigned int dim>
    class WriteToFile;
}


//should this implement the coupling::noisereduction::NoiseReduction interface?
template<unsigned int dim>
class coupling::WriteToFile : public coupling::FilterInterface<dim>{
    public:
        WriteToFile(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //covers the entire MD domain
				const std::vector<tarch::la::Vector<dim, unsigned int>> localCellIndices, //covers the entire MD domain
				std::string location):

				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices),
				_localCellIndices(localCellIndices),
		   		_location(location)
		{
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance created. Will save to: " << _location << ". Last Cell Index: " << coupling::FilterInterface<dim>::_cellIndices.back() << std::endl;
        #endif
        }

        ~WriteToFile(){
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance deconstructed." << std::endl;
        #endif
        }

     
	    void operator()();

    private:
		const std::vector<tarch::la::Vector<dim, unsigned int>> _localCellIndices;
        std::string _location;
        std::ofstream _file;
};

//include implementation of header
#include "WriteToFile.cpph"
