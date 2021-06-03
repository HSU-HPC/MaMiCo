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
 * Read-only filter that writes cell data to a specified file in .csv format.
 * @author Felix Maurer
 *
 * Information regarding the terms "mamico"/"sequence" indices, confer ../sequencing/FilterSequence.h
 *
 * Output format will be compliant to the usual MaMiCo CSV format (using ';' as separator).
 * Output order will be:
 * - current iteration
 * - sequence indexing
 * - mamico indexing
 * - scalars
 * - vectors
 *
 * The output file will either contain data of all or just the final coupling iteration. Boolean parameter 'overwrite' is used for that.
 */
template<unsigned int dim>
class coupling::WriteToFile : public coupling::FilterInterfaceReadOnly<dim>{
    public:
        WriteToFile(
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCellVector,
				const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outputCellVector,
				const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices,
				const std::array<bool, 7> filteredValues, 
				std::string location,
				bool overwrite = false,
				int oneCellOnly = -1):

				coupling::FilterInterfaceReadOnly<dim>(inputCellVector, outputCellVector, mamicoCellIndices, filteredValues, "WTF"),
				_sequenceCellIndices({}), //TODO: pass ic to this filter
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
            std::cout << "		WTF: Write to file instance created. Will save to: " << _location << ". Last Cell Index: " << coupling::FilterInterface<dim>::_cellIndices.back() << std::endl;
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
		const std::vector<tarch::la::Vector<dim, unsigned int>> _sequenceCellIndices;
        std::string _location;

		//true of only the last iteration should be in file output
		bool _overwrite;

		//-1 if all cells should be in file output, holds index of the only cell to be outputted otherwise
		int _oneCellOnly;

        std::ofstream _file;
		unsigned int _iteration; 
};

//include implementation of header
#include "WriteToFile.cpph"
