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
				bool filteredValues[7], 
				const std::vector<tarch::la::Vector<dim, unsigned int>> localCellIndices, //covers the entire MD domain
				std::string location,
				bool overwrite):

				coupling::FilterInterface<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues),
				_localCellIndices(localCellIndices),
		   		_location(location),
				_overwrite(overwrite),
				_iteration(1)
		{
			if(dim == 2 or dim == 3){
				_header += "Iteration, ";
				if(dim == 2) _header += "Global Index X, Global Index Y, ";
				else _header += "Global Index X, Global Index Y, Global Index Z, ";
				if(dim == 2) _header += "Local Index X, Local Index Y, ";
				else _header += "Local Index X, Local Index Y, Local Index Z, ";
				if(filteredValues[0]) _header += "Microscopic Mass, ";
				if(filteredValues[2]) _header += "Macroscopic Mass, "; 
				if(filteredValues[4]) _header += "Potential Energy, ";
				if(filteredValues[6]) _header += "Temperature, ";
				if(filteredValues[1]){ 			
					if(dim == 2) _header += "Microscopic Momentum X, Microscopic Momentum Y, ";
					else _header += "Microscopic Momentum X, Microscopic Momentum Y, Microscopic Momentum Z, ";
				}
				if(filteredValues[3]){
					if(dim == 2) _header += "Macroscopic Momentum X, Macroscopic Momentum Y, ";
					else _header += "Macroscopic Momentum X, Macroscopic Momentum Y, Macroscopic Momentum Z, ";
				}
				if(filteredValues[5]){
					if(dim == 2) _header += "Velocity X, Velocity Y, ";
					else _header += "Velocity X, Velocity Y, Velocity Z, ";
				}
				_header.pop_back(); //remove the last (unwanted) ' '.
				_header.pop_back(); //remove the last (unwanted) ','.
			}
			else { _header += "SINCE DIM IS NEITHER 2 NOR 3, YOU HAVE TO ADD THIS MANUALLY"; }

			if(!_overwrite){
				_file.open(location);
				_file << _header << std::endl;
				_file.close();
			}

        	#ifdef DEBUG_WRITE_TO_FILE
            std::cout << "		WTF: Write to file instance created. Will save to: " << _location << ". Last Cell Index: " << coupling::FilterInterface<dim>::_cellIndices.back() << std::endl;
			if(_overwrite) std::cout << "		It will only print output of the last iteration." << std::endl;
        	#endif
        }

        ~WriteToFile(){
        	#ifdef DEBUG_WRITE_TO_FILE
            std::cout << "		WTF: Write to file instance deconstructed." << std::endl;
        	#endif
        }

     
	    void operator()();

    private:
		const std::vector<tarch::la::Vector<dim, unsigned int>> _localCellIndices;
        std::string _location;
		bool _overwrite;

        std::ofstream _file;
		std::string _header;
		unsigned int _iteration; 
};

//include implementation of _header
#include "WriteToFile.cpph"
