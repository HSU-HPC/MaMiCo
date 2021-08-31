#pragma once

namespace coupling {
	namespace indexing {

		/*
		 * NOTE:
		 * This currently uses hardcoded values for one specific KVS test scenario.
		 * In addition to this, it is aimed to be tested sequentially, i.e. .local shall behave identical to !.local
		 *
		 * TODO:
		 * Only declarations should be made here, lower/upper boundaries must be determined at runtime using IndexingService
		 *
		 * TODO:
		 * move all of this to Service
		 */

		/*
		 * !MD TO MACRO aka MAMICO INDEXING
		 */

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

		/*
		 * MD TO MACRO
		 */

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

	}
}
