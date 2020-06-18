// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#include <string>
#include <vector>
#include <fstream>

#define DEBUG_WRITE_TO_FILE
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/filtering/FilterInterface.h"

namespace coupling {
    template<unsigned int dim>
    class WriteToFile;
}


//should this implement the coupling::noisereduction::NoiseReduction interface?
template<unsigned int dim>
class coupling::WriteToFile : public coupling::FilterInterface<dim>{
    public:
        WriteToFile(std::string location) : _location(location){
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance created. Will save to: " << _location << std::endl;
        #endif
        }

        ~WriteToFile(){
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance deconstructed." << std::endl;
        #endif
        }

     
	    void apply(
				const std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& inputCellVector,
				std::vector<coupling::datastructures::MacroscopicCell<dim>*  >& outputCellVector,
			   	const std::vector<unsigned int> cellIndices,
				const coupling::IndexConversion<dim>& indexConversion);

    private:
        std::string _location;
        std::ofstream _file;
};

//include implementation of header
#include "WriteToFile.cpph"


/*
 * Fragen:
 * Z.40 Warum wird der const Vector mit pass by reference übergeben?
 */
