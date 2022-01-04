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
	namespace filtering{
	    template<unsigned int dim>
	    class WriteToFile;
	}
}


/*
 * Read-only filter that writes cell data to a specified file in .csv format.
 *
 * Output format will be compliant to the usual MaMiCo CSV format (using ';' as separator).
 * Output order will be:
 * - current iteration
 * - scalar cell properties
 * - vector cell properties
 *
 * The output file will either contain data of all or just the final coupling iteration. Boolean parameter 'overwrite' is true in the latter case.
 * If you wish to print out just one cell and not the entire input domain, pass to "oneCellOnly" in constructor the corresponding local, md2macro, noGhost index.
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::filtering::WriteToFile : public coupling::filtering::FilterInterfaceReadOnly<dim>{
    public:
        WriteToFile(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCells,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCells,
				const std::array<bool, 7> filteredValues, 
				std::string location, //output file location
				bool overwrite = false,
				int oneCellOnly = -1):

				coupling::filtering::FilterInterfaceReadOnly<dim>(inputCells, outputCells, filteredValues, "WTF"),
		   		_location(location),
				_overwrite(overwrite),
				_oneCellOnly(oneCellOnly),
				_iteration(1)
		{	
			if(!_overwrite){
				_file.open(location);
				_file.close();
			}

        	#ifdef DEBUG_WRITE_TO_FILE
            std::cout << "		WTF: Write to file instance created. Will save to: " << _location << std::endl;
			if(_overwrite) std::cout << "		It will only print output of the last iteration." << std::endl;
			if(_oneCellOnly != -1) std::cout << "		It will only print data of cell with linear sequence domain index " << _oneCellOnly << std::endl;
        	#endif
        }

        ~WriteToFile(){
        	#ifdef DEBUG_WRITE_TO_FILE
            std::cout << "		WTF: Write to file instance deconstructed." << std::endl;
        	#endif
        }
		
	    void operator()();

    private:
        std::string _location;

		//true of only the last iteration should be in file output
		bool _overwrite;

		//-1 if all cells should be in file output, holds index of the only cell to be outputted otherwise
		int _oneCellOnly; //TODO: use CellIndex!

        std::ofstream _file;
		unsigned int _iteration; 
};

//include implementation of header
#include "WriteToFile.cpph"
