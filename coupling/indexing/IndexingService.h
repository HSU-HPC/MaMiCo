#pragma once

namespace coupling {  
	namespace indexing {

		struct IndexType{
			const bool vector = false;
			const bool local = false;
			const bool md2macro = false;
			const bool noGhost = false;

			bool constexpr operator==(const IndexType& comp) const {
				return (vector == comp.vector and local == comp.local and md2macro == comp.md2macro and noGhost == comp.noGhost);
			}
		};

		auto constexpr BaseIndexType = coupling::indexing::IndexType{true, false, false, false}; //TODO: more descriptive name
		
		template<unsigned int dim>
		class IndexingService;
	}
}

/*
 * TODO: comment
 *
 */


//Include CellIndex template class definition
#include "CellIndex.h"

//Include specialization declaration of CellIndex
#include "CellIndexDecls.h"

//Include implementation of this header. 
#include "IndexingService.cpph"

//Include non-member functions operating on indexes
#include "Operations.h"
