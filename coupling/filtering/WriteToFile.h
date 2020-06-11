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
        WriteToFile(std::string location, std::string ftype = "csv") : _location(location), _ftype(ftype){
            if(ftype != "csv"){
                std::cout << "ERROR: WriteToFile does not support file type: " << ftype << std::endl;
                exit(EXIT_FAILURE);
            }

        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance created. Will save to: " << _location << _ftype << std::endl;
        #endif
        }

        ~WriteToFile(){
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance deconstructed." << std::endl;
        #endif
        }

     
	    std::vector<coupling::datastructures::MacroscopicCell<dim>*  > apply(
				std::vector<coupling::datastructures::MacroscopicCell<dim>*  > inputCellVector,
			   	std::vector<unsigned int> cellIndices);

    private:
        std::string _location;
        std::string _ftype;
        std::ofstream _file;
};

//include implementation of header
#include "WriteToFile.cpph"


/*
 * Fragen:
 * Z.40 Warum wird der const Vector mit pass by reference übergeben?
 */
