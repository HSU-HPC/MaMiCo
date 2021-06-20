// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_ICM2M
//#define DEBUG_ICM2M_VERBOSE

#include "IndexConversion.h"
#include "interface/MacroscopicSolverInterface.h"
#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
#include <mpi.h>
#endif

namespace coupling {
	template<unsigned int dim>
	class IndexConversionMD2Macro;
}

/**
 * Wrapper class for coupling::IndexConversion for cases in which you need to know boundaries of the
 *		-- MD2Macro domain -- ,
 * that is the domain of cells that are transfered from MD to CS.
 *
 * @Author Felix Maurer
 */
template<unsigned int dim>
class coupling::IndexConversionMD2Macro {
	public:
		IndexConversionMD2Macro(
				const coupling::IndexConversion<dim>* indexConversion,
				coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface
				#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
				,const MPI_Comm comm = MPI_COMM_WORLD, //TODO: case multimd
				const int lowestRankInComm = 0 //TODO: case multimd
				#endif
			):
				_ic(indexConversion),
				_msi(macroscopicSolverInterface),
				_lowerBoundaryAllRanks(nullptr),
				_upperBoundaryAllRanks(nullptr),
				_lowerBoundaryThisRank(nullptr),
				_upperBoundaryThisRank(nullptr)
				#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
				,_comm(comm),
				_lowestRank((int)lowestRankInComm),
				_myRank(_ic != nullptr ? _ic->getThisRank() : -1)
				#endif
			{
				if(_ic == nullptr)
					throw std::runtime_error("IndexConversionMD2Macro: Constructor called with nullptr as base IndexConversion.");
				if(_msi == nullptr)
					throw std::runtime_error("IndexConversionMD2Macro: Constructor called with nullptr as MacroscopicSolverInterface.");


				#ifdef DEBUG_ICM2M
				std::cout << "ICM2M: Created new instance at location " << this << " using IC at " << _ic << " and MSI at " << _msi << std::endl;
				#endif
			}

		~IndexConversionMD2Macro() {
			delete _lowerBoundaryAllRanks;
			delete _upperBoundaryAllRanks;
			delete _lowerBoundaryThisRank; 
			delete _upperBoundaryThisRank;

			#ifdef DEBUG_ICM2M
			std::cout << "ICM2M: Deconstructed." << std::endl;
			#endif

		}

		/*
		 * This method assumes the given cells to be the entire m2m-domain (see below).
		 * It currently only has uses in cases where you do not have access to all cells of the entire domain, 
		 * which is only the case for post multi-instance filtering as of now.
		 *
		 * In this case, the "outer" region is empty.
		 */
		void initMD2MacroDomain(
				std::vector<coupling::datastructures::MacroscopicCell<dim> *> m2mDomainCells,
				std::vector<tarch::la::Vector<dim, unsigned int>> m2mGlobalCellIndices,
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outerCells,
				std::vector<tarch::la::Vector<dim, unsigned int>>& outerIndices
				);

		/*
		 * This alternative chooses a subspace of the cell (and index) input based on what will be transfered to the macro solver.
		 * This subspace is referred to as "md2Macro-domain" (or sometimes (m2m-domain).
		 * All other cells will be placed into the "outer" domain.
		 *
		 * Note that all indexing here is in terms of MaMiCo-Indexing nonetheless, i.e. m2mDomainCells don't start counting at the start of the MD-To-Macro domain.
		 */
		void initMD2MacroDomain(
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCells,
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& m2mDomainCells,
				std::vector<tarch::la::Vector<dim, unsigned int>>& m2mIndices,
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& outerCells,
				std::vector<tarch::la::Vector<dim, unsigned int>>& outerIndices
				);

		/*
		 * Writes two limiting boundary vectors to the arguments passed.
		 *
		 * For all d < dim, lowerBoundaries[d] < upperBoundaries[d].
		 * 
		 * If M2M-domain is not yet defined, this does nothing.
		 */
		void getMD2MacroDomainBoundariesAllRanks(
				tarch::la::Vector<dim, unsigned int>& lowerBoundaries,
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) const { 
			if(_lowerBoundaryAllRanks != _upperBoundaryAllRanks) {
					lowerBoundaries = *_lowerBoundaryAllRanks;
					upperBoundaries = *_upperBoundaryAllRanks;
			}
			#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
			else std::cout << "WARNING: ICM2M (" << _myRank << "): getGlobalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
			#else
			else std::cout << "WARNING: ICM2M: getGlobalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
			#endif
		}
		void getMD2MacroDomainBoundariesThisRank(
				tarch::la::Vector<dim, unsigned int>& lowerBoundaries,
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) const { 
			if(_lowerBoundaryThisRank != _upperBoundaryThisRank) {
					lowerBoundaries = *_lowerBoundaryThisRank;
					upperBoundaries = *_upperBoundaryThisRank;
			}
			#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
			else std::cout << "WARNING: ICM2M (" << _myRank << "): getLocalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
			#else
			else std::cout << "WARNING: ICM2M: getLocalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
			#endif
		}

		//assumes lower boundary to be lower than upper 
		tarch::la::Vector<dim, unsigned int> getGlobalMD2MacroDomainSize() const {
			//Since both lower and upper boundaries are inclusive, we need to add one. operator+(int) in tarch::la:Vector would be great here...
			auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
			return *_upperBoundaryAllRanks - *_lowerBoundaryAllRanks + plus_one;
		}
		tarch::la::Vector<dim, unsigned int> getLocalMD2MacroDomainSize() const {
			auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
			return *_upperBoundaryThisRank - *_lowerBoundaryThisRank + plus_one;
		}

		tarch::la::Vector<dim, double> getGlobalMD2MacroDomainOffset() const {
			tarch::la::Vector<dim, double> offset = _ic->getGlobalMDDomainOffset(); //standard MD offset

			if(_lowerBoundaryAllRanks == nullptr)
				throw std::runtime_error("ERROR: Calling while _lowerBoundaryAllRanks is uninitialized!");

			//offset of md2macro domain relative to MD domain
			for(unsigned int d = 0; d < dim; d++) {
				offset[d] += _ic->getMacroscopicCellSize()[d] * (*_lowerBoundaryAllRanks)[d]; //offset of md2macro domain relative to MD domain
			}

			//std::cout << offset << std::endl << std::endl;
			return offset;
		}

		//TODO: refactor? get -> convert

		/*
		 * Same as getGlobalVectorCellIndex but has the option to "ignore" Ghost Layer cells, i.e. return INT_MAX for all indices in Ghost Layer.
		 * E.g.: Ghost layer at x = 0. Then requesting vector index (x,y,z,..) will return (INT_MAX, y, z,...).
		 */
		tarch::la::Vector<dim,unsigned int> getGlobalVectorCellIndex(unsigned int globalCellIndex, bool noGL = true) const;
		tarch::la::Vector<dim,unsigned int> getLocalVectorCellIndex(unsigned int localCellIndex, bool noGL = true) const;

		/*
		 * In a lot of cases, you want instances of this to be indentical to a coupling::IndexConversion.
		 * This way you can access this class' underlying IC element with ease.
		 */
		const coupling::IndexConversion<dim>* getBaseIC() const { return _ic;}	

		unsigned int getLocalCellIndex(tarch::la::Vector<dim,unsigned int> localCellIndex) const{
			auto numberCells = getLocalMD2MacroDomainSize();
			unsigned int index = localCellIndex[dim-1];
			for (int d = dim-2; d >-1; d--)
				index = numberCells[d]*index + localCellIndex[d];
			return index;
		}
		
	private:

		//used during initMD2MacroDomain();
		void initGlobalMD2MacroDomainBoundaries();

		//used to determine which cells are neither withing GL or M2M domain
		const coupling::IndexConversion<dim>* _ic;
		coupling::interface::MacroscopicSolverInterface<dim>* _msi;

		/*
		 * initialised during first call of getMD2MacroDomainBoundaries
		 * these all use global indexing.
		*/
		tarch::la::Vector<dim, unsigned int>* _lowerBoundaryAllRanks;
		tarch::la::Vector<dim, unsigned int>* _upperBoundaryAllRanks;

		tarch::la::Vector<dim, unsigned int>* _lowerBoundaryThisRank;
		tarch::la::Vector<dim, unsigned int>* _upperBoundaryThisRank;

		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)	
		const MPI_Comm _comm;

		//This rank is assumed to manage cell (0,...,0) both in global and in M2M terms.
		const int _lowestRank;

		//This ICM2M instance's rank
		const int _myRank;
		#endif
};

#include "IndexConversionMD2Macro.cpph"
