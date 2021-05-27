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
//#define DEBUG_SEQ_FILTER_VERBOSE

namespace coupling{
	template<unsigned int dim>
	class SequentialFilter;
}


/*
 * WARNING: EXPERIMENTAL. TODO: output this as well during construction
 *
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
					const MPI_Comm comm = MPI_COMM_WORLD //null if run locally i.e. parallel TODO: remove default parameter, pass communicator
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

		/*
		 * These work in a very similar fashion to FilterInterface's advanced getter/setter methods.
		 */
		coupling::datastructures::MacroscopicCell<dim>* getLocalInputCellOfIndex(tarch::la::Vector<dim,unsigned int> index);
		coupling::datastructures::MacroscopicCell<dim>* getLocalOutputCellOfIndex(tarch::la::Vector<dim,unsigned int> index);

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
		 * Auxilliary functions providing an interface between low-level double buffers used by MPI and Macro Cells.
		 */
		void macroscopicCellToBuffer(std::vector<double>& buf, const coupling::datastructures::MacroscopicCell<dim>* cell);

		void bufferToMacroscopicCell(const std::vector<double>& buf, coupling::datastructures::MacroscopicCell<dim>* cell);

		void cellIndexToBuffer(std::vector<unsigned int>& buf, const tarch::la::Vector<dim, unsigned int>& index);

		void bufferToCellIndex(const std::vector<unsigned int>& buf, tarch::la::Vector<dim, unsigned int>& index);

		//The sequentialized Filter
		coupling::FilterInterface<dim>* _filter;

		//Null if run locally.
		const coupling::IndexConversionMD2Macro<dim>* _ic;	

		//MPI related stuff
		const MPI_Comm _comm;
		int _commSize;
		const int _processingRank;
		const int _myRank;

		//Globalized variants of cell and indexing data structures (i.e spanning across all cells of the global domain). Only the master rank uses these.
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells_Global; 	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Global;
		std::vector<tarch::la::Vector<dim,unsigned int>> _cellIndices_Global;

		//Buffers macro cells and indices for MPI communication
		std::vector<double> _cellbuf;
		std::vector<unsigned int> _indexbuf;

		//Used by the processing rank to remember its local domain
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells_Local;	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Local;	
		std::vector<tarch::la::Vector<dim,unsigned int>> _cellIndices_Local;
};

#include "SequentialFilter.cpph"
