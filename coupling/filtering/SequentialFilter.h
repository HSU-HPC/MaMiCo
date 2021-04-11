// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
#include "coupling/IndexConversionMD2Macro.h"
#include "coupling/CouplingMDDefinitions.h"
#include <mpi.h>


#define FILTER_SEQUENTIAL true
#define FILTER_PARALLEL false

#define DEBUG_SEQ_FILTER

namespace coupling{
	template<unsigned int dim>
	class SequentialFilter;
}


/*
 * Implementation of FilterInterface.h for filters which operate (optionally or mandatorily) in a sequential manner, i.e. process data on one master rank.
 * For such filters, operator()() will
 * 		- contribute to one dedicated processing rank: by calling contribute()
 * 		- process on master rank only and then scatter data correctly: by calling process()
 *
 * Meant to be used as a wrapper class, i.e. take a pointer to a filter object and then sequentualizes it. 
 * Used in coupling::FilterSequence<dim>::loadFiltersFromXML.
 *
 * Disclaimer: 
 * 	1. In this context 'globalized' is equivalent to 'sequentialized' and 'local' to 'parallel'.
 * 	2. Only via XML can globalized filters be added to a sequence. This implies FFF with e.g. a python function is not compatible.
 *
 * @Author Felix Maurer
 */


template <unsigned int dim>
class coupling::SequentialFilter : public coupling::FilterInterface<dim> {
	public:
		SequentialFilter(
					coupling::FilterInterface<dim>* filter,
					const coupling::IndexConversionMD2Macro<dim>* ic = nullptr, //null if run locally i.e. parallel
					const MPI_Comm comm = MPI_COMM_WORLD //null if run locally i.e. parallel TODO: in case of multimd, this goes very wrong
		);

		~SequentialFilter() {
			delete _filter;
			for(auto cell : _inputCells_Global) delete cell;
			for(auto cell : _outputCells_Global) delete cell;
	   	}
	public:
		/*
		 * Implements FilterInterface's requirement of having a ()-operand defined.
		 */
		virtual void operator()();

	private:

		/*
		 * When sequentialized, all ranks call this function.
		 */
		virtual void contribute();

		/*
		 * When sequentialized, only the processing rank calls this function. It acts as a wrapper of _filter's operator() member function.
		 */
		virtual void process(bool sequential);

		/*
		 * Auxialliary functions providing an interface between low-level double buffers used by MPI and Macro Cells.
		 */
		void macroscopicCellsToBuffer(std::vector<double>& buf, const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& cells);

		void applyBufferToMacroscopicCells(std::vector<double>& buf, const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& cells);

		//The sequentialized Filter
		coupling::FilterInterface<dim>* _filter;

		//Null if run locally.
		const coupling::IndexConversionMD2Macro<dim>* _ic;	

		//Possibility: Individual MPI communicator per globalized filter. Probably either nullptr or MPI_COMM_WORLD for now.
		const MPI_Comm _comm;
		int _commSize;
		const int _processingRank;
		const int _myRank;
		const int _cellsPerRank;

		//Globalized variants of cell and indexing data structures (i.e spanning across all cells of the global domain). Only the master rank uses these.
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells_Global; 	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Global;
		std::vector<tarch::la::Vector<dim,unsigned int>> _cellIndices_Global;

		//Buffers macro cells for MPI communication
		std::vector<double> _sendbuf;
		std::vector<double> _recvbuf;

		//only used by processing rank to keep track what to contribute/where to write output data for next filter in sequence
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells_Local;	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Local;	

		bool _firstIteration; //TODO: what is this for??
};

#include "SequentialFilter.cpph"


/*
 * TODO
 * testing
 * 	- inconsistencies in data -> always outputs data from step 0??
 * 	- segfaults if run on more than 1 rank
 */
