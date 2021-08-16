#pragma once

#include "tarch/la/Vector.h"


/* 
 * Types of cell indices:
 *	-> scalar idx VS vector<dim> idx
 *	-> MPI rank local VS global idx
 *	-> Mamico VS MD2Macro Domain idx
 *	-> NoGL idx VS total domain incl GL 
 *
 *	TODO: Proper comment
 *
 * @author Felix Maurer, Piet Jarmatz
 */


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

		// Note: this is -std=c++20
		template<unsigned int dim, IndexType idx_T = {}>
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
		CellIndex() = default;
		CellIndex(const value_T i) : _index(i){}
		
		//conversion: convert to convert_to_T
		template<coupling::indexing::IndexType convert_to_T>
		operator CellIndex<dim, convert_to_T>();

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

		// Access to primive value_T of this index
		value_T get() const { return _index; }


	private:
		static BaseIndex _lowerBoundary;
		static BaseIndex _upperBoundary;

		value_T _index;
			
};

//Include implementation of header
#include "Indexing.cpph"

//Include non-member conversion functions
#include "Conversion.h"
