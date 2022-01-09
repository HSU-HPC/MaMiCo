#pragma once


namespace coupling { 
	namespace indexing {

		/**
		 * Conversion function to convert from vector to scalar representation of CellIndex spcialisations.
		 *
		 * @tparam dim number of dimensions of the coupled simulation.
		 * @tparam idx_T index type parametrisation used by index.
		 * @param index index that shall be converted
		 * @return int representation of index
		 *
		 * @author Felix Maurer
		 */
		template<unsigned int dim, IndexTrait ... traits>
		unsigned int 
		convertToScalar(const CellIndex<dim, traits...>& index) {
			if constexpr( std::is_same_v<unsigned int, typename CellIndex<dim, traits...>::value_T> ) {
				return index.get();
			}
			else {
				//copied from deprecated coupling::IndexConversion::getCellIndex())
				
				for(unsigned d = 0; d < dim; d++)
					if(index.get()[d] < 0) throw std::runtime_error("ERROR: Indexing: Cannot convert negative vector index to scalar.");

				unsigned int i { static_cast<unsigned int>(index.get()[dim-1]) };
				for (int d = dim-2; d >-1; d--){
					i = (CellIndex<dim, traits...>::numberCellsInDomain[d])*i + static_cast<unsigned int>(index.get()[d]);
				}

				return i;
			}
		}

		/**
		 * Conversion function to convert from scalar to vector representation of CellIndex spcialisations.
		 *
		 * @tparam dim number of dimensions of the coupled simulation.
		 * @tparam idx_T index type parametrisation used by index.
		 * @param index index that shall be converted
		 * @returns vector representation of index
		 *
		 * @author Felix Maurer
		 */
		template<unsigned int dim, IndexTrait ... traits>
		tarch::la::Vector<dim, int>
		convertToVector(const CellIndex<dim, traits...>& index) {
			if constexpr( std::is_same_v<tarch::la::Vector<dim, int>, typename CellIndex<dim, traits...>::value_T>) {
				//trivial case: convert vector to vector
				return index.get();
			}
			else {
				//copied from coupling::getVectorCellIndex() 

				tarch::la::Vector<dim, int> i { 0 };
				unsigned int i_sca { index.get() };
				for (int d = dim-1; d> 0; d--){
					i[d] = static_cast<int>(i_sca / CellIndex<dim, traits...>::divisionFactor[d]);
					i_sca -= static_cast<unsigned int>(i[d])*CellIndex<dim, traits...>::divisionFactor[d];
				}
				i[0] = i_sca;

				return i;
			}
		}

	}
}
