// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterfaceReadOnly.h"

#include <string>
#include <vector>
#include <fstream>

//#define DEBUG_WRITE_TO_FILE

namespace coupling {
    template<unsigned int dim>
    class WriteToFile;
}


/*
 * TODO: Comment
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::WriteToFile : public coupling::FilterInterfaceReadOnly<dim>{
    public:
        WriteToFile(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices, //covers the entire MD domain
				const std::array<bool, 7> filteredValues, 
				const std::vector<tarch::la::Vector<dim, unsigned int>> localCellIndices, //covers the entire MD domain
				std::string location,
				bool overwrite):

				coupling::FilterInterfaceReadOnly<dim>(inputCellVector, outputCellVector, cellIndices, filteredValues, "WTF"),
				_localCellIndices(localCellIndices),
		   		_location(location),
				_overwrite(overwrite),
				_iteration(1)
		{	
			if(!_overwrite){
				_file.open(location);
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
		unsigned int _iteration; 
};

//include implementation of header
#include "WriteToFile.cpph"
