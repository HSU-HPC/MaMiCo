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
		 * NON-MD-TO-MACRO aka MAMICO INDEXING, INCL GHOST LAYER
		 */
		
		//scalar, global, !md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::divisionFactor {};

		//BaseIndex
		template<>
		BaseIndex<3> BaseIndex<3>::lowerBoundary {};
		template<>
		BaseIndex<3> BaseIndex<3>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> BaseIndex<3>::numberCellsInDomain {};
		template<>
		unsigned int BaseIndex<3>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> BaseIndex<3>::divisionFactor {};
		
		//scalar, local, !md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::local>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local>::divisionFactor {};

		//vector, local, !md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local>::divisionFactor {};


		/*
		 * MD TO MACRO, INCL GHOST LAYER
		 */

		//scalar, global, md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::md2macro>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::md2macro>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::md2macro>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro>::divisionFactor {};

		//vector, global, md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro>::divisionFactor {};

		//scalar, local, md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro>::divisionFactor {};

		//vector, local, md2macro, !noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::divisionFactor {};
		

		/*
		 * !MD TO MACRO aka MAMICO INDEXING, EXCL GHOST LAYER
		 */
		
		//scalar, global, !md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::noGhost>::divisionFactor {};

		//vector, global, !md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::noGhost>::divisionFactor {};

		//scalar, local, !md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::noGhost>::divisionFactor {};

		//vector, local, !md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::divisionFactor {};


		/*
		 * MD TO MACRO, EXCL GHOST LAYER
		 */

		//scalar, global, md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor {};

		//vector, global, md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor {};

				
		//scalar, local, md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor {};

		//vector, local, md2macro, noGL
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary {};
		template<>
		BaseIndex<3> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::numberCellsInDomain {};
		template<>
		unsigned int CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::linearNumberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::divisionFactor {};

		//declare specialisation of IndexingService
		template class IndexingService<3>;
	}

}
#else
static_assert(false, "IndexingService only available for dim=2 or dim=3.");
#endif

//impl of IndexingService
template<unsigned int dim>
void coupling::indexing::IndexingService<dim>::init(
	const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig, 
	const coupling::configurations::MaMiCoConfiguration<dim> &mamicoConfig,
	coupling::interface::MacroscopicSolverInterface<dim> * msi,
	const unsigned int rank)
{
	//init members
	_simpleMDConfig = simpleMDConfig;
	_mamicoConfig = mamicoConfig;
	_msi = msi;
	_rank = rank;

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
	CellIndex<dim>::upperBoundary = tarch::la::Vector<dim, int> { globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int> { 1 } };
	CellIndex<dim>::setDomainParameters();

	CellIndex<dim, IndexTrait::vector>::lowerBoundary = CellIndex<dim>::lowerBoundary;
	CellIndex<dim, IndexTrait::vector>::upperBoundary = CellIndex<dim>::upperBoundary;
	CellIndex<dim, IndexTrait::vector>::setDomainParameters();

	//init boundaries of all global, non-m2m, GL excluding indexing types
	CellIndex<dim, IndexTrait::noGhost>::lowerBoundary = { 1 };
	CellIndex<dim, IndexTrait::noGhost>::upperBoundary = tarch::la::Vector<dim, int> { globalNumberMacroscopicCellsInclGL - tarch::la::Vector<dim, unsigned int> { 2 } };
	CellIndex<dim, IndexTrait::noGhost>::setDomainParameters();

	CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::noGhost>::lowerBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::noGhost>::upperBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::setDomainParameters();

	//init boundaries of all global, m2m, GL excluding indexing types
	{
		CellIndex<dim> lowerBoundary { BaseIndex<dim>::lowerBoundary };
		while(not _msi->receiveMacroscopicQuantityFromMDSolver(
				tarch::la::Vector<dim, unsigned int> { CellIndex<dim, IndexTrait::vector>{ lowerBoundary }.get() })
			 ) {
			//sanity check: empty m2m domain
			if(lowerBoundary == BaseIndex<dim>::upperBoundary)
				throw std::runtime_error("IndexingService: ERROR: Empty MD-To-Macro domain!");

			//increment by one if above is too low to be in md-to-macro domain
			++lowerBoundary;
		}
		CellIndex<dim> upperBoundary { BaseIndex<dim>::upperBoundary };
		while(not _msi->receiveMacroscopicQuantityFromMDSolver(
				tarch::la::Vector<dim, unsigned int> { CellIndex<dim, IndexTrait::vector>{ upperBoundary }.get() })
			 ) {
			//sanity check: empty m2m domain 
			if(upperBoundary < lowerBoundary)
				throw std::runtime_error("IndexingService: ERROR: Empty MD-To-Macro domain!");

			//decrement by one if above is too high to be in md-to-macro domain
			--upperBoundary;
		}

		CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = lowerBoundary; 
		CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = upperBoundary; 
		CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
	}

	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

	//init boundaries of all global, m2m, GL including indexing types
	CellIndex<dim, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get() - tarch::la::Vector<dim, int> { 1 };
	CellIndex<dim, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary.get() + tarch::la::Vector<dim, int> { 1 };
	CellIndex<dim, IndexTrait::md2macro>::setDomainParameters();

	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro>::lowerBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::md2macro>::upperBoundary;
	CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::setDomainParameters();


	//handle all local indexing types
	#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //parallel scenario

		_numberProcesses = _simpleMDConfig.getMPIConfiguration().getNumberOfProcesses();

      	// determine topology offset of this rank
		const auto parallelTopologyType { _mamicoConfig.getParallelTopologyConfiguration().getParallelTopologyType() };
	  	const unsigned int scalarNumberProcesses = _numberProcesses[0] * _numberProcesses[1] * _numberProcesses[2];
      	const unsigned int parallelTopologyOffset = (_rank/scalarNumberProcesses)*scalarNumberProcesses; //copied from IndexConversion
		_parallelTopology = coupling::paralleltopology::ParallelTopologyFactory::getParallelTopology<dim>(parallelTopologyType, _numberProcesses, parallelTopologyOffset);

		std::vector<unsigned int> ranks; //used to store ranks in which certain indices occur

		//init boundaries of all local, non-m2m, GL including indexing types
		{
			BaseIndex<dim> lowerBoundary { CellIndex<dim /*global*/>::lowerBoundary }; 
			BaseIndex<dim> upperBoundary { CellIndex<dim /*global*/>::upperBoundary };
			
			//TODO: determine these two analyticaly (i.e. calculate domain bounds): you dont need to iterate over all global indices.
			while(true) {
				ranks = getRanksForGlobalIndex(lowerBoundary);
				if(std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index occurs...*/
					break;

				//sanity check: empty local domain 
				if(lowerBoundary == CellIndex<dim /*global*/>::upperBoundary) {
					using namespace std::string_literals;
					throw std::runtime_error( "IndexingService: ERROR: Empty local domain on rank "s + std::to_string(_rank) + "!"s); 
				}

				//...increment by one if above is too high to be in local domain
				++lowerBoundary; 
			}
			while(true) {
				ranks = getRanksForGlobalIndex(upperBoundary);
				if(std::find(ranks.begin(), ranks.end(), _rank) != ranks.end()) /*if _rank is found in ranks in which the tested index occurs...*/
					break;

				//sanity check: empty local domain 
				if(upperBoundary < lowerBoundary) {
					using namespace std::string_literals;
					throw std::runtime_error( "IndexingService: ERROR: Empty local domain on rank "s + std::to_string(_rank) + "!"s); 
				}

				//...decrement by one if above is too high to be in local domain
				--upperBoundary; 
			}

			CellIndex<dim, IndexTrait::local>::lowerBoundary = lowerBoundary; 
			CellIndex<dim, IndexTrait::local>::upperBoundary = upperBoundary;
			CellIndex<dim, IndexTrait::local>::setDomainParameters();
		}
		
		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::lowerBoundary = CellIndex<dim, IndexTrait::local>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::upperBoundary = CellIndex<dim, IndexTrait::local>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::setDomainParameters();

		//init boundaries of all local, non-m2m, GL excluding indexing types
		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::local>::lowerBoundary.get() + tarch::la::Vector<dim, int> { 1 };
		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::local>::upperBoundary.get() - tarch::la::Vector<dim, int> { 1 };
		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

		//init boundaries of all local, m2m, GL excluding indexing types
		{
			using LocalNoGlIndex = CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>;
			using GlobalMD2MacroNoGlIndex = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>;

			tarch::la::Vector<dim, int> lowerBoundary;
			tarch::la::Vector<dim, int> upperBoundary;
			for(unsigned int d = 0; d < dim; d++) {
				lowerBoundary[d] = std::max(LocalNoGlIndex::lowerBoundary.get()[d], GlobalMD2MacroNoGlIndex::lowerBoundary.get()[d]);
				upperBoundary[d] = std::min(LocalNoGlIndex::upperBoundary.get()[d], GlobalMD2MacroNoGlIndex::upperBoundary.get()[d]);
			}

			CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = lowerBoundary;
			CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = upperBoundary;
			CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
		}
		
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

		//init boundaries of all local, m2m, GL including indexing types
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary.get() - tarch::la::Vector<dim, int> { 1 };
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary.get() + tarch::la::Vector<dim, int> { 1 };
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

	#else //sequential scenario
		//Copy all local indexing from global
		CellIndex<dim, IndexTrait::local>::lowerBoundary = CellIndex<dim>::lowerBoundary;
		CellIndex<dim, IndexTrait::local>::upperBoundary = CellIndex<dim>::upperBoundary;
		CellIndex<dim, IndexTrait::local>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::lowerBoundary = CellIndex<dim, IndexTrait::vector>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::upperBoundary = CellIndex<dim, IndexTrait::vector>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local>::setDomainParameters();

		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>::setDomainParameters();

		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro>::lowerBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::md2macro>::upperBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::lowerBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::upperBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>::setDomainParameters();

		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();

		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::lowerBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary = CellIndex<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>::upperBoundary;
		CellIndex<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>::setDomainParameters();
	#endif
	
	#ifdef TEST_INDEXING
	//run tests
	
	//get output file for boundaries
	std::ofstream of;
	{
		using namespace std::string_literals;
		of.open("indexbounds_rank"s + std::to_string(_rank) + ".txt"s);
	}
	testing::printAllBoundaries<dim>(of);

	try {
		testing::checkTrivialConversions<dim>();
		testing::checkTrivialConversions<dim, IndexTrait::vector>();
		testing::checkTrivialConversions<dim, IndexTrait::local>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::local>();

		testing::checkTrivialConversions<dim, IndexTrait::md2macro>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::md2macro>();
		testing::checkTrivialConversions<dim, IndexTrait::local, IndexTrait::md2macro>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro>();

		testing::checkTrivialConversions<dim, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::local, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::noGhost>();

		testing::checkTrivialConversions<dim, IndexTrait::md2macro, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::md2macro, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>();
		testing::checkTrivialConversions<dim, IndexTrait::vector, IndexTrait::local, IndexTrait::md2macro, IndexTrait::noGhost>();
	}
	catch (const std::exception& e) {
		std::cout << "WARNING: IndexingService: Test indicated faulty conversion: " << e.what() << std::endl;
	}

	of.close();
	#endif
}

#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //unused in sequential scenario
/*
 * This was in large parts stolen from IndexConversion.
 */
template<unsigned int dim>
std::vector<unsigned int> 
coupling::indexing::IndexingService<dim>::getRanksForGlobalIndex(const BaseIndex<dim> &globalCellIndex) {
	std::vector<unsigned int> ranks;
	//using the old meaning of 'globalNumberMacroscopicCells' from IndexConversion
	const auto globalNumberMacroscopicCells = BaseIndex<dim>::numberCellsInDomain - tarch::la::Vector<dim, unsigned int> { 2 };

	// start and end coordinates of neighboured cells.
	tarch::la::Vector<dim,unsigned int> start(0);
	tarch::la::Vector<dim,unsigned int> end(0);
	tarch::la::Vector<dim,unsigned int> loopIndex(0);

	// determine up to 3^dim neighboured cells in the surrounding of globalCellIndex;
	// reduce this number if globalCellIndex lies on the global boundary
	for (unsigned int d = 0; d < dim; d++){
		if ((unsigned int) globalCellIndex.get()[d] > 0){start[d] = (unsigned int) globalCellIndex.get()[d]-1;}
    	end[d] = globalNumberMacroscopicCells[d] + 1;
    	if ((unsigned int) globalCellIndex.get()[d] < end[d]  ){end[d] = (unsigned int) globalCellIndex.get()[d]+1;}
  	}

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

	return ranks;
}
#endif

#if (COUPLING_MD_PARALLEL==COUPLING_MD_YES) //unused in sequential scenario
/*
 * This was in large parts stolen from IndexConversion.
 * Note that this uses the globalNumberMacroscopicCells definition excl. the ghost layer.
 */
template<unsigned int dim>
unsigned int coupling::indexing::IndexingService<dim>::getUniqueRankForMacroscopicCell(tarch::la::Vector<dim,unsigned int> globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells) const {
	//vector containing avg number of macro cells, not counting global GL. 
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
		if ( globalCellIndex[d] < averageLocalNumberMacroscopicCells[d]+1 ) {
			processCoords[d] = 0;
  		// special case: cell in last section
		} else if ( globalCellIndex[d] > averageLocalNumberMacroscopicCells[d]*(_numberProcesses[d]-1) ) {
			processCoords[d] = _numberProcesses[d]-1;
  		// all other cases
		} else {
			// remove ghost layer contribution from vector index (...-1)
			processCoords[d] = (globalCellIndex[d]-1)/averageLocalNumberMacroscopicCells[d];
		}
	}

  return _parallelTopology->getRank(processCoords); 
}
#endif
