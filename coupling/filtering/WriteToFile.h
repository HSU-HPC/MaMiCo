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
        WriteToFile(tarch::la::Vector<dim, unsigned int> domainStart,
					tarch::la::Vector<dim, unsigned int> domainEnd,
					std::string location):
					coupling::FilterInterface<dim>(domainStart, domainEnd),
		   			_location(location)
		{
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance created. Will save to: " << _location <<". Amount of Macroscopic Cells in domain (per dimension): ("<< coupling::FilterInterface<dim>::_domainSize << ")" << std::endl;
        #endif
        }

        ~WriteToFile(){
        #ifdef DEBUG_WRITE_TO_FILE
            std::cout << "WTF: Write to file instance deconstructed." << std::endl;
        #endif
        }

     
	    void apply();

    private:
        std::string _location;
        std::ofstream _file;
};

//include implementation of header
#include "WriteToFile.cpph"
