#pragma once

#include "tarch/la/Vector.h"
#include <vector>

namespace coupling {

	namespace indexing {

		/* Types of cell indices:

			-> scalar idx VS vector<dim> idx
			-> MPI rank local VS global idx
			-> Mamico VS MD2Macro Domain idx
			-> NoGL idx VS total domain incl GL 

			TODO: proper comment
		*/

		struct IndexType{
			bool vector = false;
			bool local = false;
			bool md2macro = false;
			bool noGhost = false;
		};

		// Note: this is -std=c++20
		template<unsigned int dim, IndexType idxT = {}>
		class CellIndex;
	}
}

/*
 * TODO: comment: "neutral" element of indexing type conversions
 */
auto constexpr BaseIndexType = coupling::indexing::IndexType{true, false, false, false}; //TODO: more descriptive name



/*
 * TODO: comment
 */
template<unsigned int dim, coupling::indexing::IndexType idx_T>
class coupling::indexing::CellIndex {

	/*
	 * TODO: comment
	 */
	using value_T = std::conditional_t<idx_T.vector, tarch::la::Vector<dim, unsigned int>, unsigned int>;

	using BaseIndex = CellIndex<dim, BaseIndexType>;

	public:
		//primitive constructor: entry point to indexing type system
		CellIndex(value_T i) : _index(i){}

		CellIndex(BaseIndex i) {	
			//step 2: get tarch::la::Vector representation of i and offset of this
			auto i_tarch_vec = i.get();
			auto offset = getLowerBoundary().get();
			auto ceiling = getUpperBoundary.get();

			//step 3: check if offset > i in any dimension or if i is out of bounds
			for(unsigned d = 0; d < dim; d++) {
				if(i_tarch_vec[d] < offset[d])
					throw std::runtime_error("Error: Index conversion not possible!"); //TODO: more verbose error
				if(i_tarch_vec[d] > ceiling[d])
					throw std::runtime_error("Error: Index conversion not possible!"); //TODO: more verbose error
			}

			//step 4: subtract offset and call value_T(=tarch::la::Vector) constructor
			auto result = CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(i_tarch_vec - offset); 

			//step 3: convert to scalar if neccessary
			if(!idx_T.vector)
				_index = result.swapVectorScalar().get();
			else /*idx_T.vector == true*/
				_index = result.get();
		}

		//conversion constructor: takes any other CellIndex embodying convertedType
		template<coupling::indexing::IndexType converted_T>
		CellIndex(CellIndex<dim, converted_T> i) {
			//step 1: convert i to vector if neccessary
			
			CellIndex<dim, {true, converted_T.local, converted_T.md2macro, converted_T.noGhost}> i_vec;
			if(!converted_T.vector)
				auto i_vec = i.swapVectorScalar(); 
			else
				auto i_vec = i;

			//step 2: get tarch::la::Vector representation of i and its offset
			auto i_tarch_vec = i_vec.get();
			auto converted_T_offset = CellIndex<dim, converted_T>::getLowerOffset.get();

			//step 3: check if i is out of bounds
			for(unsigned int d = 0; d < dim; d++)
				if(i_tarch_vec[d] + converted_T_offset[d] > BaseIndex::getUpperBoundary().get()[d])
					throw std::runtime_error("Error: Index conversion not possible!"); //TODO: more verbose error

			//step 4: add offset to create BaseIndex version of i
			BaseIndex base_i = BaseIndex(i_tarch_vec + converted_T_offset);

			//step 3: construct index of idx_T type from BaseIndex
			_index = CellIndex(base_i).get();
		}

		static void setBoundaries(
				BaseIndex lowerBound,
				BaseIndex upperBound
		) {
			_lowerBoundary = lowerBound;
			_upperBoundary = upperBound;
		}

		value_T get() { return _index; }


		static BaseIndex getLowerBoundary() { return _lowerBoundary; } 
		static BaseIndex getUpperBoundary() { return _upperBoundary; } 

		unsigned int operator[](unsigned int i) const { return _index[i]; } 

	protected:
		const tarch::la::Vector<dim, unsigned int> _index;

	private:
		
		//TODO: Comment, more descriptive name
		CellIndex<dim, {not idx_T.vector, idx_T.local, idx_T.md2macro, idx_T.noGhost}> 
		swapScalarVector() {
			//TODO: implement proper conversion -> copy from IndexConversion
			if(idx_T.vector) {
				return CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(42);
			}
			else /*idx_T.vector == false*/ {
				return CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(tarch::la::Vector<dim, unsigned int>(42));
			}
		}

		static BaseIndex _lowerBoundary;
		static BaseIndex _upperBoundary;

		value_T index;
			
};

//Include implementation of header
//#include "Indexing.cpph"
