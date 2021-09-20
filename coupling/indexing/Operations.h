#pragma once


namespace coupling { 
	namespace indexing {

		/**
		 * Conversion function to convert from vector to scalar representation of CellIndex spcialisations.
		 *
		 * @tparam dim number of dimensions of the coupled simulation.
		 * @tparam idx_T index type parametrisation used by index.
		 * @param index index that shall be converted
		 * @returns CellIndex with .vector=false. All other IndexType parameters remain identical.
		 *
		 * @author Felix Maurer
		 */
		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToScalar(CellIndex<dim, idx_T> index) {
			if constexpr(idx_T.vector) {
				//copied from deprecated coupling::IndexConversion::getCellIndex())
				
				unsigned int i { index.get()[dim-1] };
				for (int d = dim-2; d >-1; d--){
					i = (CellIndex<dim, idx_T>::numberCellsInDomain[d])*i + index.get()[d];
				}

				return CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}> { i };
			}
			else {
				//trivial case: convert scalar to scalar
				return index;
			}	
		}

		/**
		 * Conversion function to convert from scalar to vector representation of CellIndex spcialisations.
		 *
		 * @tparam dim number of dimensions of the coupled simulation.
		 * @tparam idx_T index type parametrisation used by index.
		 * @param index index that shall be converted
		 * @returns CellIndex with .vector=false. All other IndexType parameters remain identical.
		 *
		 * @author Felix Maurer
		 */
		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToVector(CellIndex<dim, idx_T> index) {
			if constexpr(idx_T.vector) {
				//trivial case: convert vector to vector
				return index;
			}
			else {
				//copied from coupling::getVectorCellIndex()

				tarch::la::Vector<dim,unsigned int> i { 0 };
				unsigned int i_sca { index.get() };
				for (int d = dim-1; d> 0; d--){
					i[d] = i_sca / CellIndex<dim, idx_T>::divisionFactor[d];
					i_sca -= i[d]*CellIndex<dim, idx_T>::divisionFactor[d];
				}
				i[0] = i_sca;

				return CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}> { i };
			}
		}

	}
}
