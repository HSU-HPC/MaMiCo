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
				//copied from deprecated coupling::IndexConversion::getCellIndex())
				
				//TODO: optimize: evaluate only once
				const tarch::la::Vector<dim, unsigned int> numberCells { CellIndex<dim, idx_T>::getNumberCells() };
				unsigned int i = index.get()[dim-1];
				for (int d = dim-2; d >-1; d--){
					//TODO: introduce operator[] for .vector=true CellIndex?
					i = (numberCells[d]+2)*i + index.get()[d];
				}

				return CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}> { i };
			}
			else {
				//trivial case: convert scalar to scalar
				return index;
			}	
		}


		//TODO: fixme
		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToVector(CellIndex<dim, idx_T> index) {
			if constexpr(idx_T.vector) {
				//trivial case: convert vector to vector
				return index;
			}
			else {
				//copied from coupling::getVectorCellIndex()
				
				//TODO: optimize: calculate divFactor only once
				tarch::la::Vector<dim, unsigned int> divisionFactor { 1 };
				const tarch::la::Vector<dim, unsigned int> numberCells { CellIndex<dim, idx_T>::getNumberCells() };
				for (unsigned int d = 1; d < dim; d++){
					divisionFactor[d] = divisionFactor[d-1]*(numberCells[d-1]);
				}

				std::cout << "div fac: " << divisionFactor << std::endl;

				tarch::la::Vector<dim,unsigned int> i { 0 };
				unsigned int i_sca {index.get()};
				for (int d = dim-1; d> 0; d--){
					i[d] = i_sca / divisionFactor[d];
					i_sca -= i[d]*divisionFactor[d];
				}
				i[0] = i_sca;

				return CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}> { i };
			}
		}

	}
}
