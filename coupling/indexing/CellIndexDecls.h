#pragma once

namespace coupling {
	namespace indexing {

		//BaseIndex
		template class CellIndex<3, BaseIndexType>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, BaseIndexType>::upperBoundary {0};

		//scalar, global, not m2m, noGL
		template class CellIndex<3>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3>::upperBoundary {0};

		//vector, local, not m2m, noGL
		template class CellIndex<3, {.vector=true, .local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.vector=true, .local=true}>::upperBoundary {0};

		//scalar, local, not m2m, noGL
		template class CellIndex<3, {.local=true}>;
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::lowerBoundary {0};
		template<>
		CellIndex<3, BaseIndexType> CellIndex<3, {.local=true}>::upperBoundary {0};


	}
}
