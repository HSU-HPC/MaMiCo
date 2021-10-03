// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
#include "coupling/CouplingMDDefinitions.h"
#include <mpi.h>


#define FILTER_SEQUENTIAL true
#define FILTER_PARALLEL false

//TODO: port these 
//#define DEBUG_SEQ_FILTER
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
 * 	1. In this context 'globalized' is equivalent to 'sequentialized' and 'local' to 'parallel'. TODO: fix semantics
 * 	2. Only via XML can globalized filters be added to a sequence. This implies FFF with e.g. a python function is not compatible with this feature.
 *
 * @Author Felix Maurer
 */


template <unsigned int dim>
class coupling::SequentialFilter : public coupling::FilterInterface<dim> {
	public:
		SequentialFilter(
					coupling::FilterInterface<dim>* filter, 
					const MPI_Comm comm = MPI_COMM_WORLD //TODO: remove default parameter, pass communicator
		);

		~SequentialFilter() {
			delete _filter; //TODO: possible double free?
			for(auto cell : _inputCells_Global) delete cell;
			for(auto cell : _outputCells_Global) delete cell;
	   	}
	public:
		/*
		 * Implements FilterInterface's requirement of having a ()-operand defined.
		 */
		void operator()();

		/*
		 * These work in a very similar fashion to FilterInterface's advanced getter/setter methods.
		 */
		coupling::datastructures::IndexedMacroscopicCell<dim>* getLocalInputCellOfIndex(const coupling::indexing::CellIndex<dim> index) const;
		coupling::datastructures::IndexedMacroscopicCell<dim>* getLocalOutputCellOfIndex(const coupling::indexing::CellIndex<dim> index) const;

	private:

		/*
		 * When sequentialized, all ranks call this function.
		 */
		void contribute();

		/*
		 * When sequentialized, only the processing rank calls this function. It acts as a wrapper of _filter's operator() member function.
		 */
		void process();
		
		/*
		 * Auxilliary functions providing an interface between low-level double buffers used by MPI and Macro Cells.
		 */
		void macroscopicCellToBuffer(std::vector<double>& buf, const coupling::datastructures::MacroscopicCell<dim>& cell) const;

		void bufferToMacroscopicCell(const std::vector<double>& buf, coupling::datastructures::MacroscopicCell<dim>& cell) const;

		template<coupling::indexing::IndexType idx_T>
		void cellIndexToBuffer(std::vector<unsigned int>& buf, const coupling::indexing::CellIndex<dim, idx_T>& index) const;

		template<coupling::indexing::IndexType idx_T>
		void bufferToCellIndex(const std::vector<unsigned int>& buf, coupling::indexing::CellIndex<dim, idx_T>& index) const;

		//The sequentialized Filter
		coupling::FilterInterface<dim>* _filter;

		//MPI related stuff
		const MPI_Comm _comm;
		int _commSize;
		int _processingRank;
		int _myRank;

		//Buffers macro cells and indices for MPI communication
		std::vector<double> _cellbuf;
		std::vector<unsigned int> _indexbuf;

		//Globalized variants of cell and indexing data structures (i.e spanning across all cells of the global domain). Only the master rank uses these.
		std::vector<coupling::datastructures::IndexedMacroscopicCell<dim>* > _inputCells_Global; 	
		std::vector<coupling::datastructures::IndexedMacroscopicCell<dim>* > _outputCells_Global;

		//Used by the processing rank to remember its local domain
		std::vector<coupling::datastructures::IndexedMacroscopicCell<dim>* > _inputCells_Local;	
		std::vector<coupling::datastructures::IndexedMacroscopicCell<dim>* > _outputCells_Local;	
};

#include "SequentialFilter.cpph"
