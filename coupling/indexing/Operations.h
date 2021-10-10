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
		template<unsigned int dim, bool is_vector, bool is_local, bool is_md2macro, bool is_noGhost>
		CellIndex<dim, false, is_local, is_md2macro, is_noGhost>
		convertToScalar(CellIndex<dim, is_vector, is_local, is_md2macro, is_noGhost> index) {
			if constexpr(is_vector) {
				//copied from deprecated coupling::IndexConversion::getCellIndex())
				
				unsigned int i { index.get()[dim-1] };
				for (int d = dim-2; d >-1; d--){
					i = (CellIndex<dim, is_vector, is_local, is_md2macro, is_noGhost>::numberCellsInDomain[d])*i + index.get()[d];
				}

				return CellIndex<dim, false, is_local, is_md2macro, is_noGhost> { i };
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
		template<unsigned int dim, bool is_vector, bool is_local, bool is_md2macro, bool is_noGhost>
		CellIndex<dim, true, is_local, is_md2macro, is_noGhost>
		convertToVector(CellIndex<dim, is_vector, is_local, is_md2macro, is_noGhost> index) {
			if constexpr(is_vector) {
				//trivial case: convert vector to vector
				return index;
			}
			else {
				//copied from coupling::getVectorCellIndex()

				tarch::la::Vector<dim,unsigned int> i { 0 };
				unsigned int i_sca { index.get() };
				for (int d = dim-1; d> 0; d--){
					i[d] = i_sca / CellIndex<dim, is_vector, is_local, is_md2macro, is_noGhost>::divisionFactor[d];
					i_sca -= i[d]*CellIndex<dim, is_vector, is_local, is_md2macro, is_noGhost>::divisionFactor[d];
				}
				i[0] = i_sca;

				return CellIndex<dim, true, is_local, is_md2macro, is_noGhost> { i };
			}
		}

	}
}
