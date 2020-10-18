// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once
#include "coupling/filtering/interfaces/FilterInterface.h"
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
 * 		- process and then scatter data correctly: by calling process()
 *
 * Meant to be used as a wrapper class, i.e. take a pointer to a filter object and then sequentualizes it. 
 * Used in coupling::FilterSequence<dim>::loadFiltersFromXML.
 *
 * Disclaimer: 
 * 	1. In this context 'globalized' is equivalent to 'sequentialized' and 'local' to 'parallel'.
 * 	2. Only via XML can globalized filters be added to a sequence. This implies FFF with e.g. a python function is not serializable.
 * 	3. In cases in which sequentializing is optional, a correspoding bool should be a member variable of the implementing class. All ranks should then call process(true).
 *
 * @Author Felix Maurer
 */

template <unsigned int dim>
class coupling::SequentialFilter : public coupling::FilterInterface<dim> {
	public:
		SequentialFilter(
					coupling::FilterInterface<dim>* filter,
					const coupling::IndexConversion<dim>* ic = nullptr, //null if run locally i.e. parallel
					const MPI_Comm comm = MPI_COMM_WORLD //null if run locally i.e. parallel TODO: in case of multimd, this goes very wrong
		):
			coupling::FilterInterface<dim>("SEQUENTIALIZED_FILTER"),
			_filter(filter),
			_ic(ic),
			_comm(comm),
			_processingRank(0), //TODO: case: multimd
			_cellsPerRank(_filter->getSize())
		{
			#ifdef DEBUG_SEQ_FILTER
			std::cout << "		SEQFILTER: Now initializing sequentialized Filter of type: " << _filter->getType() <<std::endl;
			#endif

			//case: sequential
			if(_ic) {
				MPI_Comm_size(comm, &_commSize);
				//allocate global cell data structures for processing rank only
				if( _processingRank == (int) _ic->getThisRank()) {
					_outputCells_Local = _filter->getOutputCells(); //filter's output cell pointers will be overwritten, but we need this information in operator()
					for(int c = 0; c < _cellsPerRank * _commSize; c++) {
						_inputCells_Global.push_back(new coupling::datastructures::MacroscopicCell<dim>());
						_outputCells_Global.push_back(new coupling::datastructures::MacroscopicCell<dim>());
					}
				}
			}
			//case: not sequential
			else {
				_commSize = 0;
				#ifdef DEBUG_SEQ_FILTER
				std::cout << "		SEQFILTER: It will not run sequentially." << std::endl;
				#endif
			}
		}

		~SequentialFilter() {
			delete _filter;
			for(auto cell : _inputCells_Global) delete cell;
			for(auto cell : _outputCells_Global) delete cell;
	   	}
	public:
		/*
		 * Implements FilterInterface's requirement of having a ()-operand defined.
		 */
		virtual void operator()(){
			if(_ic) {
				if(_processingRank == (int) _ic->getThisRank()) {
					contribute();
					process(FILTER_SEQUENTIAL); 				
				}
				else contribute();

				//Distribute output data
				MPI_Scatter(_sendbuf.data(), _cellsPerRank, MPI_DOUBLE, _recvbuf.data(), _cellsPerRank * _commSize, MPI_DOUBLE, _processingRank, _comm); 

				//Read output data from buffer
				if(_processingRank == (int) _ic->getThisRank()) applyBufferToMacroscopicCells(_recvbuf, _outputCells_Local);
				else applyBufferToMacroscopicCells(_recvbuf, _filter->getOutputCells());
			}
			else process(FILTER_PARALLEL);
		}

	private:

		/*
		 * When sequentialized, all ranks call this function.
		 */
		virtual void contribute() {
			//construct MacroscopicCell buffer from input cell vector
			macroscopicCellsToBuffer(_sendbuf, _filter->getInputCells());

			#ifdef DEBUG_SEQ_FILTER
			std::cout << "		SEQFILTER: Sending buffer of size: " << _sendbuf.size() << std::endl;
			#endif

			//allocate on processing rank
			if(_processingRank == (int) _ic->getThisRank()) _recvbuf.resize(_sendbuf.size()*_commSize);
	
			//send to processing rank 
			MPI_Gather(_sendbuf.data(), _sendbuf.size(), MPI_DOUBLE, _recvbuf.data(), _sendbuf.size() * _commSize, MPI_DOUBLE, _processingRank, _comm); 

			#ifdef DEBUG_SEQ_FILTER
			if(_processingRank == (int) _ic->getThisRank()) std::cout << "		SEQFILTER: Receiving buffer of size: " << _recvbuf.size() << std::endl;
			#endif

		}


		/*
		 * When sequentialized, only the processing rank calls this function. It acts as a wrapper of _filter's operator() member function.
		 */
		virtual void process(bool sequential) {
			//these are either local or global
			if(sequential) {

				//write all gathered cells to inputCells_Global.
				applyBufferToMacroscopicCells(_recvbuf, _inputCells_Global);

				//Get global indexing right
				//TODO for _commSize > 1
				_cellIndices_Global = _filter->getCellIndices();

				//Apply _filter
				_filter->updateCellData(_inputCells_Global, _outputCells_Global, _cellIndices_Global);
				(*_filter)();
				
				macroscopicCellsToBuffer(_sendbuf, _outputCells_Global);	
				//Now ready to scatter...
			}
			else {
				(*_filter)();
			}
		}

		void macroscopicCellsToBuffer(std::vector<double>& buf, const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& cells){
			if(cells.empty()) std::cout << "		SEQFILTER: MCTB: Warning: empty cell vector.";

			//copy cell data to buffer
			for(auto cell : cells) {
				buf.push_back(cell->getMicroscopicMass()); //std::cout << buf.back() << " ";
				for (unsigned d = 0; d < dim; d++) { buf.push_back(cell->getMicroscopicMomentum()[d]); /*std::cout << buf.back() << " ";*/ } 
				buf.push_back(cell->getMacroscopicMass()); //std::cout << buf.back() << " ";
				for (unsigned d = 0; d < dim; d++) { buf.push_back(cell->getMacroscopicMomentum()[d]); /*std::cout << buf.back() << " ";*/ } 
				buf.push_back(cell->getPotentialEnergy()); //std::cout << buf.back() << " ";
				for (unsigned d = 0; d < dim; d++) { buf.push_back(cell->getCurrentVelocity()[d]); /*std::cout << buf.back() << " ";*/ } 
				buf.push_back(cell->getTemperature()); //std::cout << buf.back() << " ";
				//std::cout << std::endl;
			}
		}	

		void applyBufferToMacroscopicCells(std::vector<double>& buf, const std::vector<coupling::datastructures::MacroscopicCell<dim> *>& cells) {

			std::cout << "		SEQFILTER: ABTMC: Buffer size: " << buf.size() << " Cell vector size: " << cells.size() << std::endl;
			if(buf.size() != cells.size() * (4+3*dim) ) throw std::runtime_error("Buffer and cell vector size must be the same!"); //TODO

			//copy buffer data to cells
			tarch::la::Vector<dim, double> mvec_buf;
			unsigned int i = 0;
			for(unsigned int c = 0; c < cells.size(); c++) {
				cells[c]->setMicroscopicMass(buf[i]); i++;
				for(unsigned int d = 0; d < dim; d++) { mvec_buf[d] = buf[i]; i++; } cells[c]->setMicroscopicMomentum(mvec_buf);
				cells[c]->setMacroscopicMass(buf[i]); i++;
				for(unsigned int d = 0; d < dim; d++) { mvec_buf[d] = buf[i]; i++; } cells[c]->setMacroscopicMomentum(mvec_buf);
				cells[c]->setPotentialEnergy(buf[i]); i++;
				for(unsigned int d = 0; d < dim; d++) { mvec_buf[d] = buf[i]; i++; } cells[c]->setCurrentVelocity(mvec_buf);
				cells[c]->setTemperature(buf[i]); i++;
			}
		}

		//The sequentialized Filter
		coupling::FilterInterface<dim>* _filter;

		//Null if run locally.
		const coupling::IndexConversion<dim>* _ic;	

		//Possibility: Individual MPI communicator per globalized filter. Probably either nullptr or MPI_COMM_WORLD for now.
		const MPI_Comm _comm;
		int _commSize;
		const int _processingRank;
		const int _cellsPerRank;

		//Globalized variants of cell and indexing data structures (i.e spanning across all cells of the global domain). Only the master rank uses these.
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCells_Global; 	
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Global;
		std::vector<tarch::la::Vector<dim,unsigned int>> _cellIndices_Global;

		//Buffers macro cells for MPI communication
		std::vector<double> _sendbuf;
		std::vector<double> _recvbuf;

		//only used by processing rank to keep track where to write output data for next filter in sequence
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _outputCells_Local;	
};


/*
 * TODO
 * testing
 */
