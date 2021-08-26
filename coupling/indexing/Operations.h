#pragma once

#include "tarch/la/Vector.h"

/*
 * TODO: comment
 *
 * @author Felix Maurer
 */

namespace coupling {
	namespace indexing {

		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToScalar(CellIndex<dim, idx_T> index) {
			if constexpr(idx_T.vector) {
				//TODO: dummy implementation, copy actual behaviour from IndexConversion
				return CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(42);
			}
			else {
				//TODO: dummy implementation, copy actual behaviour from IndexConversion
				return index;
			}	
		}



		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToVector(CellIndex<dim, idx_T> index) {
			if constexpr(idx_T.vector) {
				//TODO: dummy implementation, copy actual behaviour from IndexConversion
				return index;
			}
			else {
				//TODO: dummy implementation, copy actual behaviour from IndexConversion
				return CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(tarch::la::Vector<dim, unsigned int>(42));
			}
		}

	}
}
