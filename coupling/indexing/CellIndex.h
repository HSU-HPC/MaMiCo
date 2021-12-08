#pragma once

#include <iostream>
#include <type_traits>

#include "tarch/la/Vector.h"


namespace coupling {
	namespace indexing {

		/**
		 * Expresses type parametrisation of CellIndex specialisations.
		 *
		 * vector: implies representation as vector, false implies scalar index. \n 
		 * local: implies indexing restricted to local MD domain. \n 
		 * md2macro: implies indexing restricted to cells that are sent from MD to macro solver. \n 
		 * noGhost: implies ghost layer cells to not be included in indexing. \n 
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
			 * Checks if a given parameter pack of at least two IndexTrait entries is in the correct order.
			 * The correct order is the order in which the Traits occur in IndexTrait's declaration, i.e.:\n
			 *
			 * 		*vector < local < md2macro < noGhost*
			 *
			 * This order is naturally induced by operator> on the enum class IndexTrait.
			 *
			 * @tparam Parameter pack of at least two IndexTraits
			 * @return true iff the parameter pack is ordered correctly.
			 */
			template<IndexTrait t1, IndexTrait t2, IndexTrait ... rest>
			constexpr bool is_ordered() {
				if constexpr (sizeof...(rest) == 0) {
					return t1 < t2;
				}
				else {
					return t1 < t2 and is_ordered<t2, rest...>();
				}
			}
		}

		/**
 		* Index used to describe spatial location of a MacroscopicCell.
 		* Since various different ways of expressing this location are useful for different applications, IndexTraits are used to describe the context of this index. \n
 		*
 		* All commonly used (arithmetic) operations on MacroscopicCell indices are provided as well as seamless conversion between any two ways of expressing these indices.
 		* (cf. user-defined conversion function below)\n
 		*
 		* @tparam dim number of dimensions of the coupled simulation
 		* @tparam traits... index type parametrisation (expressed via IndexTraits) used by this specific index
		*
		* @author Felix Maurer
 		*/
		template<unsigned int dim, IndexTrait ... traits>
		class CellIndex;

		/**
		 * Base CellIndex specialisation (mainly) used for conversions.
		 */
		template<unsigned int dim>
		using BaseIndex = CellIndex<dim, IndexTrait::vector>;

		//TODO: refactor as member functions
		template<unsigned int dim, IndexTrait ... traits>
		int convertToScalar(const CellIndex<dim, traits...>&);

		template<unsigned int dim, IndexTrait ... traits>
		tarch::la::Vector<dim, int> convertToVector(const CellIndex<dim, traits...>&);
	}
}


template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
class coupling::indexing::CellIndex {
	public:

		/**
		 * Check at compile time wether traits... is in proper order, i.e.
		 *
		 * \forall i<j: traits[i] < traits[j]
		 *
		 * where '<' is the standard operator< on enum classes.
		 * This means: vector < local < md2macro < noGhost
		 * 
		 * Note that this causes duplicate IndexTraits in 'traits' to be not accepted.
		 */
		static constexpr bool checkIndexTraitOrder() {
			if constexpr(sizeof...(traits) > 1) {
				return coupling::indexing::TraitComparisons::is_ordered<traits...>();
			}
			else {
				return true;
			}
		}
		static_assert(checkIndexTraitOrder(), "Invalid order in IndexTrait parameter pack! Correct oder: IndexTrait::vector < IndexTrait::local < IndexTrait::md2macro < IndexTrait::noGhost");
		
		/**
		 * The type of this CellIndex's underlying index representation.
		 */
		using value_T = std::conditional_t<(
			coupling::indexing::TraitComparisons::is_same<coupling::indexing::IndexTrait::vector>(traits) or ...),
			tarch::la::Vector<dim, int>, 
			int
		>;

		/**
		 * Constructors
		 */
		CellIndex() = default;
		CellIndex(const CellIndex& ci) = default;
		CellIndex(const value_T& i) : _index(i){}
		
		/**
		 * Conversion function: Convert to CellIndex of same dim but different IndexType.
		 *
		 * @tparam convert_to_T IndexType parameter of the CellIndex specialisation to convert to
		 * @return CellIndex of different template parametrisation.
		 */
		template<coupling::indexing::IndexTrait ... converted_traits>
		operator CellIndex<dim, converted_traits...>() const;
	
		/**
		 * Access to primive value_T of this index.
		 * Should be used as sparingly as possible since it can lead to bugs preventable by using CellIndex instances instead of primitives.
		 *
		 * @return unsigned integer/vector representation of this index.
		 */
		value_T get() const { return _index; } //TODO: why (value_T) cast?
	
		/**
		 * Increments the index by one.
		 * Note that this does NOT increments indices in vector representation in all directions.
		 */
		CellIndex& operator++() {
			if constexpr (std::is_same_v<value_T, tarch::la::Vector<dim, int>>) {
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
			if constexpr (std::is_same_v<value_T, tarch::la::Vector<dim, int>>) {
				CellIndex<dim> scalar_base { *this };
				*this = CellIndex { --scalar_base };

			}
			else --_index;

			return *this;
		}

		/*
		 * Any two indices fulfill some relation iff the integers underlying their CellIndex<dim, {}> equivalents fulfill that relation.
		 * Note that this is different than what the static member 'contains(...)' does.
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
			numberCellsInDomain = tarch::la::Vector<dim, unsigned int> { upperBoundary.get() - lowerBoundary.get() + tarch::la::Vector<dim, int> { 1 } };
			
			linearNumberCellsInDomain = 1;
			for(unsigned int d = 0; d < dim; d++) linearNumberCellsInDomain *= numberCellsInDomain[d];
			

			tarch::la::Vector<dim, unsigned int> divFactor { 1 };
			for (unsigned int d = 1; d < dim; d++) divFactor[d] = divFactor[d-1]*(numberCellsInDomain[d-1]);
			divisionFactor = divFactor;
		}

		/**
		 * Checks if a given index is within the domain between this the lower and upper bound of this kind of indexing.
		 *
		 * @param CellIndex index index to be checked
		 * @return true iff i is within the domain of this CellIndex template specialisation.
		 */
		static bool contains(const coupling::indexing::BaseIndex<dim>& index) {
			for(unsigned int d = 0; d < dim; d++) {
				if(index.get()[d] < CellIndex::lowerBoundary.get()[d]) return false;
				if(index.get()[d] > CellIndex::upperBoundary.get()[d]) return false;
			}

			return true;
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
		 * Number of cells in this indexing's domain. Same as numberCellsInDomain, but expressed as a scalar.
		 */
		static unsigned int linearNumberCellsInDomain;

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
std::ostream& 
operator<<(std::ostream&, const coupling::indexing::CellIndex<dim, traits...>&);

/**
 * Overload arithmetic operators for CellIndex
 */
template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> 
operator+(const coupling::indexing::CellIndex<dim, traits...>&, const coupling::indexing::CellIndex<dim, traits...>&);

template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> 
operator-(const coupling::indexing::CellIndex<dim, traits...>&, const coupling::indexing::CellIndex<dim, traits...>&);

template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> 
operator+(const coupling::indexing::CellIndex<dim, traits...>&, const typename coupling::indexing::CellIndex<dim, traits...>::value_T&);

template<unsigned int dim, coupling::indexing::IndexTrait ... traits>
coupling::indexing::CellIndex<dim, traits...> 
operator-(const coupling::indexing::CellIndex<dim, traits...>&, const typename coupling::indexing::CellIndex<dim, traits...>::value_T&);

//Include implementation
#include "CellIndex.cpph"
