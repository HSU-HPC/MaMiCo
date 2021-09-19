//Include header
#include "IndexingService.h"

#include <algorithm>
#include <iterator>

/*
 * Define specialisations of CellIndex.
 * Differentiate between dim=2 and dim=3.
 *
 * Only declarations and default inits should be made here, lower/upper boundaries must be determined at runtime using IndexingService
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
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::divisionFactor {};

		//BaseIndex
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::divisionFactor {};
		
		//scalar, local, not md2macro, not noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::divisionFactor {};

		//vector, local, not md2macro, not noGL
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
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::divisionFactor {};

		//vector, global, md2macro, not noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::divisionFactor {};

		//scalar, local, md2macro, not noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::divisionFactor {};

		//vector, local, md2macro, not noGL
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
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::divisionFactor {};

		//vector, global, not md2macro, noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::divisionFactor {};

		//scalar, local, not md2macro, noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::divisionFactor {};

		//vector, local, not md2macro, noGL
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
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, global, md2macro, noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

				
		//scalar, local, md2macro, noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, local, md2macro, noGL
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

		//declare specialisation of IndexingService
		template class IndexingService<3>;
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
	coupling::interface::MacroscopicSolverInterface<dim> * msi,
	const unsigned int rank)
	: _simpleMDConfig(simpleMDConfig),
	  _mamicoConfig(mamicoConfig),
	  _msi(msi),
	  _rank(rank)
{
	//read relevant data from configs 
	const auto globalMDDomainSize { _simpleMDConfig.getDomainConfiguration().getGlobalDomainSize() }; 
	const auto macroscopicCellSize { _mamicoConfig.getMacroscopicCellConfiguration().getMacroscopicCellSize() };
	
	//calculate total number of macroscopic cells on all ranks in Base Domain
	tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells(0);
	for (unsigned int d = 0; d < dim; d++){
		globalNumberMacroscopicCells[d] = (unsigned int) floor( globalMDDomainSize[d]/macroscopicCellSize[d] + 0.5 );

		if ( fabs((globalNumberMacroscopicCells[d])*macroscopicCellSize[d] - globalMDDomainSize[d]) > 1e-13 )
			std::cout << "IndexingService: Deviation of domain size > 1e-13!" << std::endl;
	}

	//TODO: make this globalNumberMacroscopicCells and remove all usages of the old meaning (seen above)
	const auto globalNumberMacroscopicCellsInclGL { globalNumberMacroscopicCells + tarch::la::Vector<dim, unsigned int> { 2 } };

	//init boundaries of all global, non-m2m, GL including indexing types
	CellIndex<dim>::lowerBoundary = { 0 };
	CellIndex<dim>::upperBoundary = globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int> { 1 };
	CellIndex<dim>::setDomainParameters();

	CellIndex<dim, {.vector=true}>::lowerBoundary = CellIndex<dim>::lowerBoundary;
	CellIndex<dim, {.vector=true}>::upperBoundary = CellIndex<dim>::upperBoundary;
	CellIndex<dim, {.vector=true}>::setDomainParameters();

	//init boundaries of all global, non-m2m, GL excluding indexing types
	CellIndex<dim, {.noGhost=true}>::lowerBoundary = { 1 };
	CellIndex<dim, {.noGhost=true}>::upperBoundary = globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int> { 2 };
	CellIndex<dim, {.noGhost=true}>::setDomainParameters();

	CellIndex<dim, {.vector=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.noGhost=true}>::lowerBoundary;
	CellIndex<dim, {.vector=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.noGhost=true}>::upperBoundary;
	CellIndex<dim, {.vector=true, .noGhost=true}>::setDomainParameters();

	//init boundaries of all global, m2m, GL including indexing types
	CellIndex<dim> m2mGlobal_lowerBoundary { CellIndex<dim, BaseIndexType>::lowerBoundary };
	while(_msi->receiveMacroscopicQuantityFromMDSolver( CellIndex<dim, {.vector=true}>{ m2mGlobal_lowerBoundary }.get() ) == false) {
		//sanity check: empty m2m domain
		if(m2mGlobal_lowerBoundary == CellIndex<dim, BaseIndexType>::upperBoundary)
			throw std::runtime_error("IndexingService: ERROR: Empty MD-To-Macro domain!");

		//increment by one if above is too low to be in md-to-macro domain
		++m2mGlobal_lowerBoundary;
	}
	CellIndex<dim> m2mGlobal_upperBoundary { CellIndex<dim, BaseIndexType>::upperBoundary };
	while(_msi->receiveMacroscopicQuantityFromMDSolver( CellIndex<dim, {.vector=true}>{ m2mGlobal_upperBoundary }.get() ) == false) {
		//sanity check: empty m2m domain 
		if(m2mGlobal_upperBoundary < m2mGlobal_lowerBoundary)
			throw std::runtime_error("IndexingService: ERROR: Empty MD-To-Macro domain!");

		//decrement by one if above is too high to be in md-to-macro domain
		--m2mGlobal_upperBoundary;
	}

	CellIndex<dim, {.md2macro=true}>::lowerBoundary = m2mGlobal_lowerBoundary; 
	CellIndex<dim, {.md2macro=true}>::upperBoundary = m2mGlobal_upperBoundary; 
	CellIndex<dim, {.md2macro=true}>::setDomainParameters();

	CellIndex<dim, {.vector=true, .md2macro=true}>::lowerBoundary = CellIndex<dim, {.md2macro=true}>::lowerBoundary;
	CellIndex<dim, {.vector=true, .md2macro=true}>::upperBoundary = CellIndex<dim, {.md2macro=true}>::upperBoundary;
	CellIndex<dim, {.vector=true, .md2macro=true}>::setDomainParameters();

	//init boundaries of all global, m2m, GL excluding indexing types
	//note that m2m overrules GL by definition, i.e. .noGhost has no effect if .md2macro == true
	CellIndex<dim, {.md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.md2macro=true}>::lowerBoundary;
	CellIndex<dim, {.md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.md2macro=true}>::upperBoundary;
	CellIndex<dim, {.md2macro=true, .noGhost=true}>::setDomainParameters();

	CellIndex<dim, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.md2macro=true}>::lowerBoundary;
	CellIndex<dim, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.md2macro=true}>::upperBoundary;
	CellIndex<dim, {.vector=true, .md2macro=true, .noGhost=true}>::setDomainParameters();


	//handle all local indexing types
	#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //parallel scenario

		_numberProcesses = _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses(); //TODO read this properly

      	// determine topology offset of this rank
		const auto parallelTopologyType { _mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType() };
	  	const unsigned int scalarNumberProcesses = _numberProcesses[0] * _numberProcesses[1] * _numberProcesses[2];
      	const unsigned int parallelTopologyOffset = (_rank/scalarNumberProcesses)*scalarNumberProcesses; //copied from IndexConversion
		_parallelTopology = coupling::paralleltopology::ParallelTopologyFactory::getParallelTopology<dim>(parallelTopologyType, _numberProcesses, parallelTopologyOffset);

		std::vector<unsigned int> ranks; //used to store ranks in which certain indices occur

		//init boundaries of all local, non-m2m, GL including indexing types
		CellIndex<dim, BaseIndexType> local_lowerBoundary { CellIndex<dim /*global*/>::lowerBoundary }; //used to test which indices are within local bounds
		while(true) {
			ranks = getRanksForGlobalIndex(local_lowerBoundary, globalNumberMacroscopicCells);
			if(std::ranges::find(ranks, _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index occurs...*/
				break;

			//sanity check: empty local domain 
			if(local_lowerBoundary == CellIndex<dim /*global*/>::upperBoundary) {
				using namespace std::string_literals;
				throw std::runtime_error( "IndexingService: ERROR: Empty local domain on rank "s + std::to_string(_rank) + "!"s); 
			}

			//...increment by one if above is too high to be in md-to-macro domain
			++local_lowerBoundary; 
		}
		CellIndex<dim, BaseIndexType> local_upperBoundary { CellIndex<dim /*global*/>::upperBoundary };
		while(true) {
			ranks = getRanksForGlobalIndex(local_upperBoundary, globalNumberMacroscopicCells);
			if(std::ranges::find(ranks, _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index occurs...*/
				break;

			//sanity check: empty local domain 
			if(local_upperBoundary < local_lowerBoundary) {
				using namespace std::string_literals;
				throw std::runtime_error( "IndexingService: ERROR: Empty local domain on rank "s + std::to_string(_rank) + "!"s); 
			}

			//...decrement by one if above is too high to be in md-to-macro domain
			--local_upperBoundary; 
		}

		CellIndex<dim, {.local=true}>::lowerBoundary = local_lowerBoundary; 
		CellIndex<dim, {.local=true}>::upperBoundary = local_upperBoundary;
		CellIndex<dim, {.local=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true}>::lowerBoundary = CellIndex<dim, {.vector=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true}>::upperBoundary = CellIndex<dim, {.vector=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true}>::setDomainParameters();

		//init boundaries of all local, non-m2m, GL excluding indexing types
		CellIndex<dim, {.local=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.local=true}>::lowerBoundary.get() + tarch::la::Vector<dim, unsigned int> { 1 };
		CellIndex<dim, {.local=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.local=true}>::upperBoundary.get() - tarch::la::Vector<dim, unsigned int> { 1 };
		CellIndex<dim, {.local=true, .noGhost=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.local=true, .noGhost=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.local=true, .noGhost=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::setDomainParameters();

		//init boundaries of all local, m2m, GL including indexing types
		CellIndex<dim, {.local=true}> m2mLocal_lowerBoundary { CellIndex<dim, {.local=true}>::lowerBoundary };
		while(_msi->receiveMacroscopicQuantityFromMDSolver( CellIndex<dim, BaseIndexType>{ m2mLocal_lowerBoundary }.get() ) == false) {
			//sanity check: empty m2m domain 
			if(m2mLocal_lowerBoundary == CellIndex<dim, {.local=true}>::upperBoundary) {
				std::cout << "IndexingService: WARNING: Empty local MD-To-Macro domain on rank " << _rank << "!" << std::endl;
				break; 
			}

			//increment by one if above is too high to be in md-to-macro domain
			++m2mLocal_lowerBoundary; 
		}
		CellIndex<dim, {.local=true}> m2mLocal_upperBoundary { CellIndex<dim, {.local=true}>::upperBoundary };
		while(_msi->receiveMacroscopicQuantityFromMDSolver( CellIndex<dim, BaseIndexType>{ m2mLocal_upperBoundary }.get() ) == false) {
			//sanity check: empty m2m domain
			if(m2mLocal_upperBoundary < m2mLocal_lowerBoundary) {
				std::cout << "IndexingService: WARNING: Empty local MD-To-Macro domain on rank " << _rank << "!" << std::endl;
				break;
			}

			//decrement by one if above is too high to be in md-to-macro domain
			--m2mLocal_upperBoundary; 
		}

		CellIndex<dim, {.local=true, .md2macro=true}>::lowerBoundary = m2mLocal_lowerBoundary;
		CellIndex<dim, {.local=true, .md2macro=true}>::upperBoundary = m2mLocal_upperBoundary;
		CellIndex<dim, {.local=true, .md2macro=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary = CellIndex<dim, {.local=true, .md2macro=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::upperBoundary = CellIndex<dim, {.local=true, .md2macro=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::setDomainParameters();

		//init boundaries of all local, m2m, GL excluding indexing types
		//note that m2m overrules GL by definition, i.e. .noGhost has no effect if .md2macro == true
		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.local=true, .md2macro=true}>::lowerBoundary;
		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.local=true, .md2macro=true}>::upperBoundary;
		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();

	#else //sequential scenario
		//Copy all local indexing from global
		CellIndex<dim, {.local=true}>::lowerBoundary = CellIndex<dim>::lowerBoundary;
		CellIndex<dim, {.local=true}>::upperBoundary = CellIndex<dim>::upperBoundary;
		CellIndex<dim, {.local=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true}>::lowerBoundary = CellIndex<dim, {.vector=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true}>::upperBoundary = CellIndex<dim, {.vector=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true}>::setDomainParameters();

		CellIndex<dim, {.local=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.noGhost=true}>::lowerBoundary;
		CellIndex<dim, {.local=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.noGhost=true}>::upperBoundary;
		CellIndex<dim, {.local=true, .noGhost=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.vector=true, .noGhost=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.vector=true, .noGhost=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .noGhost=true}>::setDomainParameters();

		CellIndex<dim, {.local=true, .md2macro=true}>::lowerBoundary = CellIndex<dim, {.md2macro=true}>::lowerBoundary;
		CellIndex<dim, {.local=true, .md2macro=true}>::upperBoundary = CellIndex<dim, {.md2macro=true}>::upperBoundary;
		CellIndex<dim, {.local=true, .md2macro=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary = CellIndex<dim, {.vector=true, .md2macro=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::upperBoundary = CellIndex<dim, {.vector=true, .md2macro=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true}>::setDomainParameters();

		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.md2macro=true, .noGhost=true}>::lowerBoundary;
		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.md2macro=true, .noGhost=true}>::upperBoundary;
		CellIndex<dim, {.local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();

		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary = CellIndex<dim, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary = CellIndex<dim, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary;
		CellIndex<dim, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::setDomainParameters();
	#endif
}

/*
 * This was in large parts stolen from IndexConversion.
 */
template<unsigned int dim>
std::vector<unsigned int> 
coupling::indexing::IndexingService<dim>::getRanksForGlobalIndex(const CellIndex<dim, BaseIndexType> &globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells) {
	std::vector<unsigned int> ranks;

	// start and end coordinates of neighboured cells.
	tarch::la::Vector<dim,unsigned int> start(0);
	tarch::la::Vector<dim,unsigned int> end(0);
	tarch::la::Vector<dim,unsigned int> loopIndex(0);

	// determine up to 3^dim neighboured cells in the surrounding of globalCellIndex;
	// reduce this number if globalCellIndex lies on the global boundary
	for (unsigned int d = 0; d < dim; d++){
		if (globalCellIndex.get()[d] > 0){start[d] = globalCellIndex.get()[d]-1;}
    	end[d] = globalNumberMacroscopicCells[d] + 1;
    	if (globalCellIndex.get()[d] < end[d]  ){end[d] = globalCellIndex.get()[d]+1;}
  	}

	/*
	 * TODO: refactor
	 * FM: I dont really get what is going on in {-- .. --}.
	 */
	// {--
	// loop over neighbouring regions
	for (loopIndex[2] = start[2]; loopIndex[2] <= end[2]; loopIndex[2]++){
		for (loopIndex[1] = start[1]; loopIndex[1] <= end[1]; loopIndex[1]++){
			for (loopIndex[0] = start[0]; loopIndex[0] <= end[0]; loopIndex[0]++){

				// determine the global cell index of this particular grid cell
				tarch::la::Vector<dim,unsigned int> thisGlobalCellIndex(0);
				for (unsigned int d = 0; d <dim; d++){thisGlobalCellIndex[d] = loopIndex[d];}

				// determine the unique rank for this cell
				const unsigned int rank = getUniqueRankForMacroscopicCell(thisGlobalCellIndex, globalNumberMacroscopicCells); 

				// add this rank to the vector with all ranks if we did not add this one before
				bool isContained = false;
				const unsigned int thisSize = (unsigned int) ranks.size();
				for (unsigned int i = 0; i < thisSize; i++){
					if (ranks[i] == rank){isContained=true; break;}
				}
				if (!isContained){ranks.push_back(rank);}
			}
		}
	}
	// --}

	return ranks;
}

//TODO: inline everything below here
/*
 * This was in large parts stolen from IndexConversion.
 * Note that this uses the globalNumberMacroscopicCells definition excl. the ghost layer.
 */
template<unsigned int dim>
unsigned int coupling::indexing::IndexingService<dim>::getUniqueRankForMacroscopicCell(tarch::la::Vector<dim,unsigned int> globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells) const {

	//vector containing avg number of macro cells, not counting global GL. //TODO: make this a member?
	tarch::la::Vector<dim, unsigned int> averageLocalNumberMacroscopicCells { 0 };
	for (unsigned int d = 0; d < dim; d++) {
		if(globalCellIndex[d] >= globalNumberMacroscopicCells[d]+2) { //greater or equal to the total global number incl GL (+2)
			using namespace std::string_literals;
			throw std::runtime_error("IndexingService: getUniqueRankForMacroscopicCell(): Global cell index greater than global size in dim "s + std::to_string(d));
		}
			
		averageLocalNumberMacroscopicCells[d] = globalNumberMacroscopicCells[d]/_numberProcesses[d];
	}

	tarch::la::Vector<dim,unsigned int> processCoords(0);
	for (unsigned int d = 0; d < dim; d++){
		// special case: cell in first section
		if ( globalCellIndex[d] < averageLocalNumberMacroscopicCells[d]+1 ) [[unlikely]] {
			processCoords[d] = 0;
  		// special case: cell in last section
		} else if ( globalCellIndex[d] > averageLocalNumberMacroscopicCells[d]*(_numberProcesses[d]-1) ) [[unlikely]] {
			processCoords[d] = _numberProcesses[d]-1;
  		// all other cases
		} else [[likely]]{
			// remove ghost layer contribution from vector index (...-1)
			processCoords[d] = (globalCellIndex[d]-1)/averageLocalNumberMacroscopicCells[d];
		}
	}

  return _parallelTopology->getRank(processCoords); 
}
