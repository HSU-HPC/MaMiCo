#pragma once

#include <iostream>
#include <type_traits>

#include "tarch/la/Vector.h"


namespace coupling {
	namespace indexing {

		/**
		 * Stores type parametrisation of CellIndex specialisation.
		 *
		 * .vector: True implies representation as vector, false implies scalar index. \n 
		 * .local: True implies indexing restricted to local MD domain. \n 
		 * .md2macro: True implies indexing restricted to cells that are sent from MD to macro solver. \n 
		 * .noGhost: True implies ghost layer cells to not be included in indexing. \n 
		 *
		 * @author Felix Maurer, Piet Jarmatz
		 */
		struct IndexType{
			const bool vector = false;
			const bool local = false;
			const bool md2macro = false;
			const bool noGhost = false;

			/**
			 * Pointwise equality operator for IndexType instances.
			 *
			 * @param comp IndexType to compare *this to.
			 * @returns true iff all four boolean parameters match between the two IndexType instances
			 */
			bool constexpr operator==(const IndexType& comp) const {
				return (vector == comp.vector and local == comp.local and md2macro == comp.md2macro and noGhost == comp.noGhost);
			}
		};

		/**
		 * This specific instance of IndexType is used to enable streamlined conversions and computations on CellIndex objects.
		 */
		auto constexpr BaseIndexType = coupling::indexing::IndexType{true, false, false, false}; 

		// Note: this is -std=c++20
		template<unsigned int dim, IndexType idx_T = {}>
		class CellIndex;
	}
}


/**
 * Index used to describe spatial location of a MacroscopicCell.
 * Since various different ways of expressing this location are useful for different applications, IndexType is used to describe the context of this index.
 *
 * All commonly used (arithmetic) operations on MacroscopicCell indices are provided as well as seamless conversion between any two ways of expressing these indices.
 * (cf. user-defined conversion function below)
 *
 * @tparam dim number of dimensions of the coupled simulation
 * @tparam idx_T index type parametrisation used by this specific index
 */
template<unsigned int dim, coupling::indexing::IndexType idx_T>
class coupling::indexing::CellIndex {
	public:

		/**
		 * The type of this CellIndex's underlying index representation.
		 */
		using value_T = std::conditional_t<idx_T.vector, tarch::la::Vector<dim, unsigned int>, unsigned int>;
		/**
		 * CellIndex specialisation using BaseIndexType for this cell's number of dimensions.
		 */
		using BaseIndex = CellIndex<dim, BaseIndexType>;

		//primitive constructors
		CellIndex() = default;
		CellIndex(const value_T i) : _index(i){}
		CellIndex(const CellIndex& ci) : _index(ci.get()){}
		
		/**
		 * Conversion function: Convert to CellIndex of same dim but different IndexType.
		 *
		 * @tparam convert_to_T IndexType parameter of the CellIndex specialisation to convert to
		 * @returns CellIndex of different template parametrisation.
		 */
		template<coupling::indexing::IndexType convert_to_T>
		operator CellIndex<dim, convert_to_T>() const;
	
		/**
		 * Access to primive value_T of this index.
		 * Should be used as sparingly as possible since it can lead to bugs preventable by using CellIndex instances instead of primitives.
		 *
		 * @returns unsigned integer/vector representation of this index.
		 */
		value_T get() const { return (value_T) _index; }

		//friend functions: overload arithmetic operators
		friend CellIndex operator+<>(const CellIndex &i1, const CellIndex &i2);
		friend CellIndex operator-<>(const CellIndex &i1, const CellIndex &i2);

		/**
		 * Increments the index by one.
		 * Note that this does NOT increments indices in vector representation in all directions.
		 */
		CellIndex& operator++() {
			if constexpr (idx_T.vector) {
				CellIndex<dim> scalar_base { *this };
				*this = CellIndex { ++scalar_base };
			}
			else ++_index;

			return *this;
		}
		/**
		 * Decrements the index by one.
		 * Note that this does NOT decrements indices in vector representation in all directions.
		 */
		CellIndex& operator--() {
			if constexpr (idx_T.vector) {
				CellIndex<dim> scalar_base { *this };
				*this = CellIndex { --scalar_base };

			}
			else --_index;

			return *this;
		}

		/*
		 * Any two indices fulfill some relation iff the unsigned integers underlying their CellIndex<dim, {}> equivalents fulfill that relation.
		 *
		 * @param CellIndex index to compare this index to
		 */
		bool operator==(const CellIndex &) const = default;
		bool operator!=(const CellIndex &) const = default;
		bool operator<(const CellIndex &i) const { return ( convertToScalar<dim, idx_T>(*this).get() < convertToScalar<dim, idx_T>(i).get() ); };
		bool operator<=(const CellIndex &i) const { return ( convertToScalar<dim, idx_T>(*this).get() <= convertToScalar<dim, idx_T>(i).get() ); };
		bool operator>(const CellIndex &i) const { return ( convertToScalar<dim, idx_T>(*this).get() > convertToScalar<dim, idx_T>(i).get() ); };
		bool operator>=(const CellIndex &i) const { return ( convertToScalar<dim, idx_T>(*this).get() >= convertToScalar<dim, idx_T>(i).get() ); };


		/**
		 * Initialises all static members dependant only on upperBoundary and lowerBoundary
		 */
		static void setDomainParameters() {
			numberCellsInDomain = upperBoundary.get() - lowerBoundary.get() + tarch::la::Vector<dim, unsigned int> { 1 };

			tarch::la::Vector<dim, unsigned int> divFactor { 1 };
			for (unsigned int d = 1; d < dim; d++) divFactor[d] = divFactor[d-1]*(numberCellsInDomain[d-1]);
			divisionFactor = divFactor;
		}

		/**
		 * Defines where this type of cell index starts counting. 
		 * Read inclusively, e.g.: lowerBoundary = {1,2,3} means {1,2,3} is the first index contained in this type of cell index' domain.
		 */
		static BaseIndex lowerBoundary;
		/**
		 * Defines where this type of cell index stops counting. 
		 * Read inclusively, e.g.: upperBoundary = {4,5,6} means {4,5,6} is the last index contained in this type of cell index' domain.
		 */
		static BaseIndex upperBoundary;

		/**
		 * Number of cells in this indexing's domain. Because the above declared boundaries are inclusive, this is never 0 in any direction.
		 * Initialised in setDomainParameters().
		 */
		static tarch::la::Vector<dim, unsigned int> numberCellsInDomain;
		/**
		 * Used in scalar -> vector indexing conversion functions
		 * Initialised in setDomainParameters().
		 */
		static tarch::la::Vector<dim, unsigned int> divisionFactor;

	private:
		value_T _index;

};

/**
 * Overload operator<< for CellIndex
 */
template<unsigned int dim, coupling::indexing::IndexType idx_T>
std::ostream& operator<<(std::ostream& os, const coupling::indexing::CellIndex<dim, idx_T>& i);


//Include implementation
#include "CellIndex.cpph"
