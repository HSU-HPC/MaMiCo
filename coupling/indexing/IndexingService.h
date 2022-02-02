#pragma once

#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
// parallel topologies
#include "coupling/paralleltopology/ParallelTopology.h"
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

//Include CellIndex template class definition
#include "CellIndex.h"

namespace coupling {
	namespace indexing {

		template<unsigned int dim>
		class IndexingService;

		template<unsigned int dim>
		std::vector<unsigned int> 
		getRanksForGlobalIndex(const BaseIndex<dim> &globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells);

	}
}

//Include non-member functions operating on indexes
#include "Operations.h"

//enable/disable tests
//#define TEST_INDEXING

#ifdef TEST_INDEXING
//Inlcude index tests
#include "Testing.h"
#endif

/**
 * Singleton service class initialising lower and upper boundaries of all possible CellIndex specialisations.
 *
 * @tparam dim number of dimensions of the coupled simulation
 * @param simpleMDConfig config object of SimpleMD instance used in coupling
 * @param mamicoConfig config object containg general information of coupling process
 * @param msi pointer to interface of coupled macroscopic solver
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::indexing::IndexingService{
	public:
		static IndexingService& getInstance() {
			static IndexingService singleton {};
			return singleton;
		}

		void init(	const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig,
					const coupling::configurations::MaMiCoConfiguration<dim> &mamicoConfig,
					coupling::interface::MacroscopicSolverInterface<dim> *msi,
					const unsigned int rank);

		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //parallel scenario
		/**
		 * Determines all ranks that contain a certain global BaseIndex.
		 * Ripped from deprecated IndexConversion.
		 *
		 * @param globalCellIndex index to be looked up
		 * @param globalNumberMacroscopicCells global number of cells in BaseIndex domain EXCLUDING global ghost layer cells.
		 * @returns vector of all cells which contain the index
		 */
		std::vector<unsigned int> getRanksForGlobalIndex(const BaseIndex<dim> &globalCellIndex) const;
		#endif

		unsigned int getRank() const {return _rank;}

	private:
		#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //parallel scenario
		/**
		 * Helper function used by getRanksForGlobalIndex().
		 */
		//TODO inline in getRanksForGlobalIndex()
		unsigned int getUniqueRankForMacroscopicCell(tarch::la::Vector<dim,unsigned int> globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells) const;

		/*const*/ tarch::la::Vector<dim, unsigned int> _numberProcesses; //TODO: make const
		const coupling::paralleltopology::ParallelTopology<dim> *_parallelTopology;
		#endif

		simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
		coupling::configurations::MaMiCoConfiguration<dim> _mamicoConfig;
		coupling::interface::MacroscopicSolverInterface<dim> *_msi; 
		unsigned int _rank;
};
