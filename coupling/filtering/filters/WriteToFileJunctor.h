// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_WTF_JUNCTION

#include "coupling/filtering/interfaces/JunctorInterface.h"
#include "coupling/filtering/filters/WriteToFile.h"

namespace coupling {
  template<unsigned int dim>
  class WriteToFileJunctor;
}

/** 
 * Combines two WritToFile objects into one FilterJunctor.
 * You can use this e.g. to output both primary and secondary cell data in an AsymmetricalFilterJunction.
 *
 * TODO: allow aribitrary amounts of wtfs to be joint?
 *
 * @author Felix Maurer
 * 
 */
template<unsigned int dim>
class coupling::WriteToFileJunctor : public coupling::JunctorInterface<dim,2,2> {
	public:
		WriteToFileJunctor(//first cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector1,
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector2,
			const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices1,
			const std::vector<tarch::la::Vector<dim, unsigned int>> sequenceCellIndices1,
			//second cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector1,
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector2,
			const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices2,
			//"global" parameters for both WriteToFile instances
			const std::array<bool, 7> filteredValues,

			//WriteToFile-specific parameters. [0] is for the first WriteToFile instance and [1] for the second one respectively.
			std::array<std::string,2> location,
			std::array<bool,2> overwrite = { false },
			std::array<int,2> oneCellOnly = { -1 }):
		  
			coupling::JunctorInterface<dim,2,2>( 
				{ inputCellVector1, inputCellVector2 },
		   		{ outputCellVector1, outputCellVector2 }, 
				{},//We dont really ever need either of the two cellIndices later. This is rather unclean, but JunctorInterface asks for a set of indices.
				filteredValues,
		   		"WTF-J")
		{
			_writeToFile1 = new coupling::WriteToFile<dim>(	inputCellVector1,
															outputCellVector1,
															mamicoCellIndices1,
															sequenceCellIndices1,
															filteredValues,
															location[0],
															overwrite[0],
															oneCellOnly[0]);
			_writeToFile2 = new coupling::WriteToFile<dim>(	inputCellVector2,
															outputCellVector2,
															cellIndices2,
															{},
															filteredValues,
															location[1],
															overwrite[1],
															oneCellOnly[1]);
		}
      
		~WriteToFileJunctor() {
			delete _writeToFile1;
			delete _writeToFile2;
			#ifdef DEBUG_WTF_JUNCTION
				std::cout << "    WTF-J: Destroyed WriteToFileJunctor instance." << std::endl;
			#endif
		}

    void operator()() {
		#ifdef DEBUG_WTF_JUNCTION
				std::cout << "    WTF-J: Now calling operator() on both WriteToFile instances." << std::endl;
		#endif
		(*_writeToFile1)();
		(*_writeToFile2)();
	}

  private:
	//Both WTF instances are created during construction and deleted in destruction
	coupling::WriteToFile<dim>* _writeToFile1;
	coupling::WriteToFile<dim>* _writeToFile2;
};

