//Include header
#include "IndexingService.h"

/*
 * Define specialisations of CellIndex.
 * Differentiate between dim=2 and dim=3.
 *
 * TODO:
 * Only declarations should be made here, lower/upper boundaries must be determined at runtime using IndexingService
 */

#ifdef MDDim2
//TODO

#elif MDDim3
namespace coupling {
	namespace indexing {

		/*
		 * Declare specialisations of CellIndex.
		 * Define their static members.
		 */


		/*
		 * !MD TO MACRO aka MAMICO INDEXING, INCL GHOST LAYER
		 */
		
		//scalar, global, not md2macro, not noGL
		template class CellIndex<3>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::divisionFactor {};

		//BaseIndex
		template class CellIndex<3, BaseIndexType>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::divisionFactor {};
		
		//scalar, local, not md2macro, not noGL
		template class CellIndex<3, {.local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::divisionFactor {};

		//vector, local, not md2macro, not noGL
		template class CellIndex<3, {.vector=true, .local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true}>::divisionFactor {};


		/*
		 * MD TO MACRO, INCL GHOST LAYER
		 */

		//scalar, global, md2macro, not noGL
		template class CellIndex<3, {.md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::divisionFactor {};

		//vector, global, md2macro, not noGL
		template class CellIndex<3, {.vector=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::divisionFactor {};

		//scalar, local, md2macro, not noGL
		template class CellIndex<3, {.local=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::divisionFactor {};

		//vector, local, md2macro, not noGL
		template class CellIndex<3, {.vector=true, .local=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::divisionFactor {};
		

		/*
		 * !MD TO MACRO aka MAMICO INDEXING, EXCL GHOST LAYER
		 */
		
		//scalar, global, not md2macro, noGL
		template class CellIndex<3, {.noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::divisionFactor {};

		//vector, global, not md2macro, noGL
		template class CellIndex<3, {.vector=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::divisionFactor {};

		//scalar, local, not md2macro, noGL
		template class CellIndex<3, {.local=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::divisionFactor {};

		//vector, local, not md2macro, noGL
		template class CellIndex<3, {.vector=true, .local=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::divisionFactor {};


		/*
		 * MD TO MACRO, EXCL GHOST LAYER
		 */

		//scalar, global, md2macro, noGL
		template class CellIndex<3, {.md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, global, md2macro, noGL
		template class CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

				
		//scalar, local, md2macro, noGL
		template class CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, local, md2macro, noGL
		template class CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

		//declare specialisation of IndexingService
		template<>
		class IndexingService<3>;
	}

}
#else
static_assert(false, "IndexingService only available for dim=2 or dim=3.");

#endif

//impl of IndexingService
template<unsigned int dim>
coupling::indexing::IndexingService<dim>::IndexingService(
	const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig, 
	const coupling::configurations::MaMiCoConfiguration<dim> &mamicoConfig,
	const coupling::interface::MacroscopicSolverInterface<dim> *msi)
	: _simpleMDConfig(simpleMDConfig),
	  _mamicoConfig(mamicoConfig),
	  _msi(msi)
		//TODO: align to the left
		{
			//read relevant data from configs 
			const auto globalMDDomainSize { _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize() };
			const auto macroscopicCellSize { _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize() };
			
			//calculate total number of macroscopic cells on all ranks in Base Domain
			tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells(0);
			for (unsigned int d = 0; d < dim; d++){
				globalNumberMacroscopicCells[d] = (unsigned int) floor( globalMDDomainSize[d]/macroscopicCellSize[d] + 0.5 );

				if ( fabs(globalNumberMacroscopicCells[d]*macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13 )
					std::cout << "IndexingService: Deviation of domain size > 1e-13!" << std::endl;
			}

			//init boundaries of all global, non-m2m, GL including indexing types
			CellIndex<3>::lowerBoundary = 0;
			CellIndex<3>::upperBoundary = globalNumberMacroscopicCells - tarch::la::Vector<dim, unsigned int> {1};
			CellIndex<3>::setDomainParameters();

			CellIndex<3, {.vector=true}>::lowerBoundary = CellIndex<3>::lowerBoundary;
			CellIndex<3, {.vector=true}>::upperBoundary = CellIndex<3>::upperBoundary;
			CellIndex<3, {.vector=true}>::setDomainParameters();

			//init boundaries of all global, non-m2m, GL excluding indexing types
			CellIndex<3, {.noGhost=true}>::lowerBoundary = 1;
			CellIndex<3, {.noGhost=true}>::upperBoundary = globalNumberMacroscopicCells - tarch::la::Vector<dim, unsigned int> {2};
			CellIndex<3, {.noGhost=true}>::setDomainParameters();

			CellIndex<3, {.vector=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.noGhost=true}>::lowerBoundary;
			CellIndex<3, {.vector=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.noGhost=true}>::upperBoundary;
			CellIndex<3, {.vector=true, .noGhost=true}>::setDomainParameters();

			//init boundaries of all global, m2m, GL including indexing types
			auto m2mGlobal_lowerBoundary { CellIndex<3, BaseIndexType>::lowerBoundary };
			while(_msi->receiveMacroscopicQuantityFromSolver( m2mGlobal_lowerBoundary.get() ) == false)
				//increment by one (in all dims) if above is too low to be in md-to-macro domain
				m2mGlobal_lowerBoundary += tarch::la::Vector<dim, unsigned int> { 1 };
			auto m2mGlobal_upperBoundary { CellIndex<3, BaseIndexType>::upperBoundary };
			while(_msi->receiveMacroscopicQuantityFromSolver( m2mGlobal_upperBoundary.get() ) == false)
				//decrement by one (in all dims) if above is too high to be in md-to-macro domain
				m2mGlobal_upperBoundary -= tarch::la::Vector<dim, unsigned int> { 1 };

			CellIndex<3, {.md2macro=true}>::lowerBoundary = m2mGlobal_lowerBoundary; 
			CellIndex<3, {.md2macro=true}>::upperBoundary = m2mGlobal_upperBoundary; 
			CellIndex<3, {.md2macro=true}>::setDomainParameters();

			CellIndex<3, {.vector=true, .md2macro=true}>::lowerBoundary = CellIndex<3, {.md2macro=true}>::lowerBoundary;
			CellIndex<3, {.vector=true, .md2macro=true}>::upperBoundary = CellIndex<3, {.md2macro=true}>::upperBoundary;
			CellIndex<3, {.vector=true, .md2macro=true}>::setDomainParameters();

			//init boundaries of all global, m2m, GL excluding indexing types
			//note that m2m overrules GL by definition, i.e. .noGhost has no effect if .md2macro == true
			CellIndex<3, {.md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.md2macro=true}>::lowerBoundary;
			CellIndex<3, {.md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.md2macro=true}>::upperBoundary;
			CellIndex<3, {.md2macro=true, .noGhost=true}>::setDomainParameters();

			CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.md2macro=true}>::lowerBoundary;
			CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.md2macro=true}>::upperBoundary;
			CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::setDomainParameters();


			//handle all local indexing types
			#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //parallel scenario
				static_assert(false, "ERROR: Indexing system does not yet work for parallel scenarios!");

				//read more data from configs
				const auto innerOverlap { _mamicoConfig.getMomentumInsertionConfiguration().getInnerOverlap() };
				/* TODO
				const auto globalMDDomainOffset { };
				const auto numberProcesses { };
				const unsigned int rank { };

				const auto parallelTopologyType { };
				unsigned int topologyOffset { };*/

			#else //sequential scenario
				//Copy all local indexing from global
				CellIndex<3, {.local=true}>::lowerBoundary = CellIndex<3>::lowerBoundary;
				CellIndex<3, {.local=true}>::upperBoundary = CellIndex<3>::upperBoundary;
				CellIndex<3, {.local=true}>::setDomainParameters();

				CellIndex<3, {.vector=true, .local=true}>::lowerBoundary = CellIndex<3, {.vector=true}>::lowerBoundary;
				CellIndex<3, {.vector=true, .local=true}>::upperBoundary = CellIndex<3, {.vector=true}>::upperBoundary;
				CellIndex<3, {.vector=true, .local=true}>::setDomainParameters();

				CellIndex<3, {.local=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.noGhost=true}>::lowerBoundary;
				CellIndex<3, {.local=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.noGhost=true}>::upperBoundary;
				CellIndex<3, {.local=true, .noGhost=true}>::setDomainParameters();

				CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.vector=true, .noGhost=true}>::lowerBoundary;
				CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.vector=true, .noGhost=true}>::upperBoundary;
				CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::setDomainParameters();

				CellIndex<3, {.local=true, .md2macro=true}>::lowerBoundary = CellIndex<3, {.md2macro=true}>::lowerBoundary;
				CellIndex<3, {.local=true, .md2macro=true}>::upperBoundary = CellIndex<3, {.md2macro=true}>::upperBoundary;
				CellIndex<3, {.local=true, .md2macro=true}>::setDomainParameters();

				CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary = CellIndex<3, {.vector=true, .md2macro=true}>::lowerBoundary;
				CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::upperBoundary = CellIndex<3, {.vector=true, .md2macro=true}>::upperBoundary;
				CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::setDomainParameters();

				CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.md2macro=true, .noGhost=true}>::lowerBoundary;
				CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.md2macro=true, .noGhost=true}>::upperBoundary;
				CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();

				CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary;
				CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary;
				CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();
			#endif
		}

