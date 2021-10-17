// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#include "coupling/filtering/interfaces/JunctorInterface.h"
#include "coupling/filtering/filters/Datastructures.h"

#define AVG_JUNCTOR_DEBUG

namespace coupling {
	template<unsigned int dim>
	class AveragingJunctor;
}

/** 
 * Simple Junctor (cf. interfaces/FilterJunctor.h) averaging cell data of two input sets.
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::AveragingJunctor : public coupling::JunctorInterface<dim,2,1> {
public:
	using coupling::JunctorInterface<dim,2,1>::_inputCellVectors;
	using coupling::JunctorInterface<dim,2,1>::_outputCellVectors;
	using coupling::FilterInterface<dim>::_scalarGetters;
	using coupling::FilterInterface<dim>::_vectorGetters;
	using coupling::FilterInterface<dim>::_scalarSetters;
	using coupling::FilterInterface<dim>::_vectorSetters;

	AveragingJunctor(
		  const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector1,
		  const std::vector<coupling::datastructures::MacroscopicCell<dim> *> inputCellVector2,
		  const std::vector<coupling::datastructures::MacroscopicCell<dim> *> outputCellVector,
		  const std::vector<tarch::la::Vector<dim, unsigned int>> cellIndices,
		  const std::array<bool, 7> filteredValues
	):
		coupling::JunctorInterface<dim,2,1>( 
				{ inputCellVector1 , inputCellVector2 },
				{ outputCellVector }, 
				cellIndices, 
				filteredValues,
				"Average")
	{}

	virtual ~AveragingJunctor(){
		#ifdef AVG_JUNCTOR_DEBUG
		std::cout << "    AVG: Averaging Junctor instance destructed." << std::endl;
		#endif
	}

	void operator()() {
		/*
		 * Fill output set by averaging cell data from the two input cell sets.
		 */
		for(unsigned int ci = 0; ci < _outputCellVectors[0].size(); ci++) {
			//Iterate over all scalar values that are filtered.
			for(unsigned int si = 0; si < _scalarSetters.size(); si++) {
				(_outputCellVectors[0]/*first (and only) output partition*/[ci]->*(_scalarSetters[si])) //call setter on output cell
					(  (_inputCellVectors[0]/*first input partition*/[ci]->*(_scalarGetters[si]))() * 0.5 //call getter on cell from input set 1
					 + (_inputCellVectors[1]/*second input partition*/[ci]->*(_scalarGetters[si]))() * 0.5); //call getter on cell from input set 2
			}

			//Iterate over all multidimensional values that are filtered.
			for(unsigned int vi = 0; vi < _vectorSetters.size(); vi++) {
				(_outputCellVectors[0]/*first (and only) output partition*/[ci]->*(_vectorSetters[vi])) //call setter on output cell
					(  (_inputCellVectors[0]/*first input partition*/[ci]->*(_vectorGetters[vi]))() * 0.5 //call getter on cell from input set 1
					 + (_inputCellVectors[1]/*second input partition*/[ci]->*(_vectorGetters[vi]))() * 0.5); //call getter on cell from input set 2
			}
		}
	}

};
