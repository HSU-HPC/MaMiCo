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
		//primitive constructors
		CellIndex() : _index() {}
		CellIndex(const value_T i) : _index(i){}
		
		//conversion: convert to convert_to_T
		template<coupling::indexing::IndexType convert_to_T>
		operator CellIndex<dim, convert_to_T>();

		value_T get() const { return _index; }

		/*
		 * This must be initialised for each index type individually.
		 * Together with the type of value_t, these two boundaries define a CellType.
		 */
		static void setBoundaries(
				BaseIndex lowerBound,
				BaseIndex upperBound
		) {
			_lowerBoundary = lowerBound;
			_upperBoundary = upperBound;
		}

		static BaseIndex getLowerBoundary() { return _lowerBoundary; } 
		static BaseIndex getUpperBoundary() { return _upperBoundary; } 

		unsigned int operator[](unsigned int i) const { return _index[i]; } 
		
		//TODO: Comment
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


	private:
		static BaseIndex _lowerBoundary;
		static BaseIndex _upperBoundary;

		value_T _index;
			
};

//Include implementation of header
#include "Indexing.cpph"
