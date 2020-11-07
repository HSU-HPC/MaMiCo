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
 * TODO: Extensive comments
 * @Author Felix Maurer
 */
template<unsigned int dim>
//TODO: Is this really supposed to inherit from IC?
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
				_lowestRank(lowestRankInComm)
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
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) { 
			if(_globalLowerBoundaries != _globalUpperBoundaries) {
					lowerBoundaries = *_globalLowerBoundaries;
					upperBoundaries = *_globalUpperBoundaries;
			}
			else {/*TODO: Warning. This is only the case if both are NULL.*/}
		}
		void getLocalMD2MacroDomainBoundaries(
				tarch::la::Vector<dim, unsigned int>& lowerBoundaries,
				tarch::la::Vector<dim, unsigned int>& upperBoundaries) { 
			if(_localLowerBoundaries != _localUpperBoundaries) {
					lowerBoundaries = *_localLowerBoundaries;
					upperBoundaries = *_localUpperBoundaries;
			}
			else {/*TODO: Warning. See above.*/}
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
		const coupling::IndexConversion<dim>* operator()() const { return _ic;}	

		
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
		const unsigned int _lowestRank;
};

#include "IndexConversionMD2Macro.cpph"
