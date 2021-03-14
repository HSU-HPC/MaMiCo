// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_ICM2M

#include "IndexConversion.h"
#include "interface/MacroscopicSolverInterface.h"
#include <mpi.h>

namespace coupling {
	template<unsigned int dim>
	class IndexConversionMD2Macro;
}

/**
 * Wrapper class for coupling::IndexConversion for cases in which you need to know boundaries of the
 *		-- MD2Macro domain -- ,
 * that is the domain of cells that are transfered from MD to CS.
 * To do so, it makes use of coupling::services::MacroscopicCellService and MPI communcation between processes.
 *
 * @Author Felix Maurer
 */
template<unsigned int dim>
class coupling::IndexConversionMD2Macro {
	public:
		IndexConversionMD2Macro(
				const coupling::IndexConversion<dim>* indexConversion,
				coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface,
				const MPI_Comm comm = MPI_COMM_WORLD, //TODO: case multimd
				const int lowestRankInComm = 0 //TODO: case multimd
			):
				_ic(indexConversion),
				_msi(macroscopicSolverInterface),
				_globalLowerBoundaries(nullptr),
				_globalUpperBoundaries(nullptr),
				_localLowerBoundaries(nullptr),
				_localUpperBoundaries(nullptr),
				_comm(comm),
				_lowestRank((int)lowestRankInComm),
				_myRank((int)_ic->getThisRank())
			{
				#ifdef DEBUG_ICM2M
				std::cout << "ICM2M: Created new instance at location " << this << " using IC at " << _ic << " and MSI at " << _msi << std::endl;
				#endif
			}

		~IndexConversionMD2Macro() {
			delete _globalLowerBoundaries;
			delete _globalUpperBoundaries;
			delete _localLowerBoundaries; 
			delete _localUpperBoundaries;

			#ifdef DEBUG_ICM2M
			std::cout << "ICM2M: Deconstructed." << std::endl;
			#endif

		}

		/*
		 * Chooses a subspace of the cell (and index) input based on what will be transfered to the macro solver.
		 * This subspace is referred to as "md2Macro-domain" (or sometimes (m2m-domain).
		 */
		void initMD2MacroDomain(
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& inputCells,
				std::vector<coupling::datastructures::MacroscopicCell<dim> *>& m2mDomainCells,
				std::vector<tarch::la::Vector<dim, unsigned int>>& m2mIndices
				);

		/*
		 * Writes two limiting boundary vectors to the arguments passed.
		 *
		 * For all d < dim, lowerBoundaries[d] < upperBoundaries[d].
		 * 
		 * If M2M-domain is not yet defined, this does nothing.
		 */
		void getGlobalMD2MacroDomainBoundaries(
				tarch::la::Vector<dim, unsigned int>& lowerBoundaries,
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) const { 
			if(_globalLowerBoundaries != _globalUpperBoundaries) {
					lowerBoundaries = *_globalLowerBoundaries;
					upperBoundaries = *_globalUpperBoundaries;
			}
			else std::cout << "WARNING: ICM2M (" << _myRank << "): getGlobalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
		}
		void getLocalMD2MacroDomainBoundaries(
				tarch::la::Vector<dim, unsigned int>& lowerBoundaries,
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) const { 
			if(_localLowerBoundaries != _localUpperBoundaries) {
					lowerBoundaries = *_localLowerBoundaries;
					upperBoundaries = *_localUpperBoundaries;
			}
			else std::cout << "WARNING: ICM2M (" << _myRank << "): getLocalMD2MacroDomainBoundaries while domain boundaries are unitialized!" << std::endl; 
		}

		//assumes lower boundary to be lower than upper 
		tarch::la::Vector<dim, unsigned int> getGlobalMD2MacroDomainSize() const {
			//Since both lower and upper boundaries are inclusive, we need to add one. operator+(int) in tarch::la:Vector would be great here...
			auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
			return *_globalUpperBoundaries - *_globalLowerBoundaries + plus_one;
		}
		tarch::la::Vector<dim, unsigned int> getLocalMD2MacroDomainSize() const {
			auto plus_one = tarch::la::Vector<dim, unsigned int>(1);
			return *_localUpperBoundaries - *_localLowerBoundaries + plus_one;
		}

		tarch::la::Vector<dim, double> getGlobalMD2MacroDomainOffset() const {
			auto offset = _ic->getGlobalMDDomainOffset(); //standard MD offset

			//offset of md2macro domain relative to MD domain
			for(unsigned int d = 0; d < dim; d++) 
				offset[d] += _ic->getMacroscopicCellSize()[d] * (*_globalLowerBoundaries)[d]; //offset of md2macro domain relative to MD domain

			//std::cout << offset << std::endl << std::endl;
			return offset;
		}


		//TODO: move both of these to regular IC?
		/*
		 * Same as getGlobalVectorCellIndex but sets all ghost layer indices to INT_MAX
		 * E.g.: Ghost layer at x = 0. Then requesting vector index (x,y,z,..) will return (INT_MAX, y, z,...).
		 */
		tarch::la::Vector<dim,unsigned int> getGlobalVectorCellIndex_NoGL(unsigned int globalCellIndex) const;
		tarch::la::Vector<dim,unsigned int> getLocalVectorCellIndex_NoGL(unsigned int localCellIndex) const;

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
		//used to determine which cells are neither withing GL or M2M domain
		const coupling::IndexConversion<dim>* _ic;
		coupling::interface::MacroscopicSolverInterface<dim>* _msi;

		//initialised during first call of getMD2MacroDomainBoundaries
		tarch::la::Vector<dim, unsigned int>* _globalLowerBoundaries;
		tarch::la::Vector<dim, unsigned int>* _globalUpperBoundaries;

		tarch::la::Vector<dim, unsigned int>* _localLowerBoundaries;
		tarch::la::Vector<dim, unsigned int>* _localUpperBoundaries;

		const MPI_Comm _comm;

		//This rank is assumed to manage cell (0,...,0) both in global and in M2M terms.
		const int _lowestRank;

		//This ICM2M instance's rank
		const int _myRank;
};

#include "IndexConversionMD2Macro.cpph"
