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
		convertToScalar(CellIndex<dim, idx_T> index);

		template<unsigned int dim, IndexType idx_T>
		CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
		convertToVector(CellIndex<dim, idx_T> index);
	}

}

using namespace coupling::indexing;

/*
 * Pattern matching convertToScalar:
 * Set nullptr as default parameter in order to match template signature of forward decls found above.
 */
template<unsigned int dim, IndexType idx_T, 
	//Pattern: idx_T is a vector
	std::enable_if_t<idx_T.vector>* = nullptr >
CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
convertToScalar(CellIndex<dim, idx_T> index) {
	//TODO: dummy implementation, copy actual behaviour from IndexConversion
	return CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(42);
}

template<unsigned int dim, IndexType idx_T, 
	//Pattern: idx_T is no vector, i.e. a scalar
	std::enable_if_t<not idx_T.vector>* = nullptr >
CellIndex<dim, {false, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
convertToScalar(CellIndex<dim, idx_T> index) {
	return index;
}

/*
 * Pattern matching convertToVector:
 */
template<unsigned int dim, IndexType idx_T, 
	//Pattern: idx_T is a vector
	std::enable_if_t<idx_T.vector>* = nullptr >
CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
convertToVector(CellIndex<dim, idx_T> index) {
	return index;
}

template<unsigned int dim, IndexType idx_T, 
	//Pattern: idx_T is no vector, i.e. a scalar
	std::enable_if_t<not idx_T.vector>* = nullptr >
CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>
convertToVector(CellIndex<dim, idx_T> index) {
	//TODO: dummy implementation, copy actual behaviour from IndexConversion
	return CellIndex<dim, {true, idx_T.local, idx_T.md2macro, idx_T.noGhost}>(tarch::la::Vector<dim, unsigned int>(42));
}

//Overloading operator<< for CellIndex
template<unsigned int dim, coupling::indexing::IndexType idx_T>
std::ostream& operator<<(std::ostream& os, const CellIndex<dim, idx_T>& i) {
	os << i.get();
	return os;
}

//Overloading arithmetic operators for CellIndex using friend functions
template<unsigned int dim, IndexType idx_T>
CellIndex<dim, idx_T> operator+(
		const CellIndex<dim, idx_T>& i1,
		const CellIndex<dim, idx_T>& i2
) { return CellIndex<dim, idx_T>(i1._index + i2._index); }

template<unsigned int dim, coupling::indexing::IndexType idx_T>
coupling::indexing::CellIndex<dim, idx_T> operator-(
		const CellIndex<dim, idx_T>& i1,
		const CellIndex<dim, idx_T>& i2
) { return CellIndex<dim, idx_T>(i1._index - i2._index); }
