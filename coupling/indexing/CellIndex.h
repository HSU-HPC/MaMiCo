#pragma once

#include <iostream>
#include <type_traits>
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
	
		// Note: this is -std=c++20
		template<unsigned int dim, IndexType idx_T = {}>
		class CellIndex;
	}
}


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
		CellIndex(const CellIndex& ci) : _index(ci.get()){}
		
		//conversion: convert to convert_to_T
		template<coupling::indexing::IndexType convert_to_T>
		operator CellIndex<dim, convert_to_T>() const;
	
		//access to primive value_T of this index
		value_T get() const { return (value_T) _index; }
		//explicit operator value_T() const { return _index; } TODO: enable this?

		//friend functions: overload arithmetic operators
		friend CellIndex operator+<>(const CellIndex &i1, const CellIndex &i2);
		friend CellIndex operator-<>(const CellIndex &i1, const CellIndex &i2);

		//TODO: move init of boundaries here. decide where this should be called. constexpr -> Decls?
		static void setDomainParameters() {
			numberCellsInDomain = upperBoundary.get() - lowerBoundary.get() + tarch::la::Vector<dim, unsigned int> { 1 };

			tarch::la::Vector<dim, unsigned int> divFactor { 1 };
			for (unsigned int d = 1; d < dim; d++) divFactor[d] = divFactor[d-1]*(numberCellsInDomain[d-1]);
			divisionFactor = divFactor;
		}

		/*
		 * Note: both are inclusive
		 */
		static BaseIndex lowerBoundary;
		static BaseIndex upperBoundary;

		//Number of cells in this indexing's domain. Because the above declared boundaries are inclusive, this is never 0 in any direction.
		static tarch::la::Vector<dim, unsigned int> numberCellsInDomain;
		//Used in scalar -> vector indexing functions
		static tarch::la::Vector<dim, unsigned int> divisionFactor;

	private:
		const value_T _index;

};

//overload operator<<
template<unsigned int dim, coupling::indexing::IndexType idx_T>
std::ostream& operator<<(std::ostream& os, const coupling::indexing::CellIndex<dim, idx_T>& i);



//Include implementation of header
#include "CellIndex.cpph"
