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
		 * NOTE:
		 * This currently uses hardcoded values for one specific KVS test scenario.
		 * In addition to this, it is aimed to be tested sequentially, i.e. .local shall behave identical to !.local
		 */


		/*
		 * !MD TO MACRO aka MAMICO INDEXING, INCL GHOST LAYER
		 */
		
		//scalar, global, not md2macro, not noGL
		template class CellIndex<3>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::upperBoundary {13};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3>::divisionFactor {};

		//BaseIndex
		template class CellIndex<3, BaseIndexType>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::upperBoundary {13};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, BaseIndexType>::divisionFactor {};
		
		//scalar, local, not md2macro, not noGL
		template class CellIndex<3, {.local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::upperBoundary {13};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true}>::divisionFactor {};

		//vector, local, not md2macro, not noGL
		template class CellIndex<3, {.vector=true, .local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::upperBoundary {13};
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
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true}>::divisionFactor {};

		//vector, global, md2macro, not noGL
		template class CellIndex<3, {.vector=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::numberCellsInDomain {0};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true}>::divisionFactor {0};

		//scalar, local, md2macro, not noGL
		template class CellIndex<3, {.local=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true}>::divisionFactor {};

		//vector, local, md2macro, not noGL
		template class CellIndex<3, {.vector=true, .local=true, .md2macro=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true}>::upperBoundary {9};
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
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::lowerBoundary {1};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.noGhost=true}>::upperBoundary {12};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.noGhost=true}>::divisionFactor {};

		//vector, global, not md2macro, noGL
		template class CellIndex<3, {.vector=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::lowerBoundary {1};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .noGhost=true}>::upperBoundary {12};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .noGhost=true}>::divisionFactor {};

		//scalar, local, not md2macro, noGL
		template class CellIndex<3, {.local=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::lowerBoundary {1};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .noGhost=true}>::upperBoundary {12};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .noGhost=true}>::divisionFactor {};

		//vector, local, not md2macro, noGL
		template class CellIndex<3, {.vector=true, .local=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::lowerBoundary {1};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .noGhost=true}>::upperBoundary {12};
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
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.md2macro=true, .noGhost=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, global, md2macro, noGL
		template class CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {0};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .md2macro=true, .noGhost=true}>::divisionFactor {0};

				
		//scalar, local, md2macro, noGL
		template class CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

		//vector, local, md2macro, noGL
		template class CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::lowerBoundary {4};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::upperBoundary {9};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::numberCellsInDomain {};
		template<>
		tarch::la::Vector<3, unsigned int> CellIndex<3, {.vector=true, .local=true, .md2macro=true, .noGhost=true}>::divisionFactor {};

	}
}
#else
static_assert(false, "IndexingService only available for dim=2 or dim=3.");

#endif
