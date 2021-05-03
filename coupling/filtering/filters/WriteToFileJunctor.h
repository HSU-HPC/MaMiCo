// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#define DEBUG_WTF_JUNCTION

#include "coupling/filtering/interfaces/AsymmetricalJunctorInterface.h"
#include "coupling/filtering/filters/WriteToFile.h"

namespace coupling {
  template<unsigned int dim>
  class WriteToFileJunctor;
}

/** 
 * Combines two WriteToFile objects into one FilterJunctor.
 * You can use this e.g. to output both primary and secondary cell data in an AsymmetricalFilterJunction.
 * TODO: more explanation
 *
 * @author Felix Maurer
 * 
 */
template<unsigned int dim>
class coupling::WriteToFileJunctor : public coupling::AsymmetricalJunctorInterface<dim> {
	public:
		WriteToFileJunctor(
			//first cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector1,
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector1,
			const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices1,
			const std::vector<tarch::la::Vector<dim, unsigned int>> sequenceCellIndices1,

			//second cell data set
			const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector2,
			//no output cells
			const std::vector<tarch::la::Vector<dim, unsigned int>> mamicoCellIndices2,
			//no sequence indices
	
			//"global" parameters for both WriteToFile instances
			const std::array<bool, 7> filteredValues,

			//WriteToFile-specific parameters. [0] is for the first WriteToFile instance and [1] for the second one respectively.
			std::array<std::string,2> location,
			std::array<bool,2> overwrite = { false },
			std::array<int,2> oneCellOnly = { -1 }):
		  
			coupling::AsymmetricalJunctorInterface<dim>( 
				inputCellVector1,
				outputCellVector1,
				mamicoCellIndices1,
				
				inputCellVector2,
				mamicoCellIndices2,
				filteredValues,
		   		"WTF-J")
		{
			//write to file instance covering first cell data set
			coupling::AsymmetricalJunctorInterface<dim>::_filter1 
				= new coupling::WriteToFile<dim>(	inputCellVector1,
													outputCellVector1,
													mamicoCellIndices1,
													sequenceCellIndices1,
													filteredValues,
													location[0],
													overwrite[0],
													oneCellOnly[0]);

			//write to file instance covering second cell data set
			coupling::AsymmetricalJunctorInterface<dim>::_filter2
				= new coupling::WriteToFile<dim>(	inputCellVector2,
													{}, //no output
													mamicoCellIndices2,
													{}, //no sequence indexing
													filteredValues,
													location[1],
													overwrite[1],
													oneCellOnly[1]);
		}
      
		~WriteToFileJunctor() {
			#ifdef DEBUG_WTF_JUNCTION
				std::cout << "    WTF-J: Destroyed WriteToFileJunctor instance." << std::endl;
			#endif
		}
};

