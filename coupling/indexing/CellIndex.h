#pragma once

#include <iostream>
#include <type_traits>

#include "tarch/la/Vector.h"


/**
 * TODO:
 * update comments for cpp17 port
 */

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
		 * @author Felix Maurer
		 */
		enum class IndexTrait {vector, local, md2macro, noGhost};
		namespace TraitComparisons {

			/**
			 * Returns true iff the template and runtime argument match. Curried operator== for above enum class.
			 */
			template<IndexTrait t1>
			constexpr bool is_same(const IndexTrait& t2) { return t1 == t2; }

			/**
			 * Returns true iff the template pack contains t.
			 */
			//template<IndexTrait ... Ts>
			//constexpr bool contains(const IndexTrait t) { return (is_same<t>(Ts) or ...; }
		}

		//TODO comment
		template<unsigned int dim, IndexTrait ... traits>
		class CellIndex;

		/**
		 * Base CellIndex specialisation (mainly) used for conversions.
		 */
		template<unsigned int dim>
		using BaseIndex = CellIndex<dim, IndexTrait::vector>;

		template<unsigned int dim, IndexTrait ... traits>
		unsigned int convertToScalar(const CellIndex<dim, traits...>&);

		template<unsigned int dim, IndexTrait ... traits>
		tarch::la::Vector<dim, unsigned int> convertToVector(const CellIndex<dim, traits...>&);
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
 *
 *
 *TODO: problem: different ordering in Ts creates different types -> enable_if usage?
 */
template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
class coupling::indexing::CellIndex {
	public:

		/**
		 * The type of this CellIndex's underlying index representation.
		 */
		using value_T = std::conditional_t<(
			coupling::indexing::TraitComparisons::is_same<coupling::indexing::IndexTrait::vector>(traits) or ...),
			tarch::la::Vector<dim, unsigned int>, 
			unsigned int
		>;

		//primitive constructors
		CellIndex() = default;
		CellIndex(const value_T& i) : _index(i){}
		CellIndex(const CellIndex& ci) : _index(ci.get()){}
		
		/**
		 * Conversion function: Convert to CellIndex of same dim but different IndexType.
		 *
		 * @tparam convert_to_T IndexType parameter of the CellIndex specialisation to convert to
		 * @returns CellIndex of different template parametrisation.
		 */
		template<coupling::indexing::IndexTrait ... converted_traits>
		operator CellIndex<dim, converted_traits...>() const;
	
		/**
		 * Access to primive value_T of this index.
		 * Should be used as sparingly as possible since it can lead to bugs preventable by using CellIndex instances instead of primitives.
		 *
		 * @returns unsigned integer/vector representation of this index.
		 */
		value_T get() const { return (value_T) _index; }
	
		/**
		 * Increments the index by one.
		 * Note that this does NOT increments indices in vector representation in all directions.
		 */
		CellIndex& operator++() {
			if constexpr (std::is_same_v<value_T, tarch::la::Vector<dim, unsigned int>>) {
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
			if constexpr (std::is_same_v<value_T, tarch::la::Vector<dim, unsigned int>>) {
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
		bool operator==(const CellIndex& i) const { return _index == i.get(); }
		bool operator!=(const CellIndex& i) const { return not (i == *this); } 
		bool operator<(const CellIndex& i) const { return convertToScalar<dim, traits...>(*this) < convertToScalar<dim, traits...>(i); };
		bool operator<=(const CellIndex& i) const { return convertToScalar<dim, traits...>(*this) <= convertToScalar<dim, traits...>(i); };
		bool operator>(const CellIndex& i) const { return convertToScalar<dim, traits...>(*this) > convertToScalar<dim, traits...>(i); };
		bool operator>=(const CellIndex& i) const { return convertToScalar<dim, traits...>(*this) >= convertToScalar<dim, traits...>(i); };


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
		static BaseIndex<dim> lowerBoundary;
		/**
		 * Defines where this type of cell index stops counting. 
		 * Read inclusively, e.g.: upperBoundary = {4,5,6} means {4,5,6} is the last index contained in this type of cell index' domain.
		 */
		static BaseIndex<dim> upperBoundary;

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
template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
std::ostream& operator<<(std::ostream&, const coupling::indexing::CellIndex<dim, traits...>&);

/**
 * Overload arithmetic operators for CellIndex
 */
template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> operator+(const coupling::indexing::CellIndex<dim, traits...>&, const coupling::indexing::CellIndex<dim, traits...>&);

template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> operator-(const coupling::indexing::CellIndex<dim, traits...>&, const coupling::indexing::CellIndex<dim, traits...>&);


//Include implementation
#include "CellIndex.cpph"
