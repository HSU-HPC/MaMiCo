// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#pragma once

#define DEBUG_FILTER_JUNCTION_ASYM

#include "coupling/filtering/sequencing/FilterSequence.h"

//INCLUDE ALL RELEVANT JUNCTOR HEADERS HERE
//TODO: implement (asymmetrical) junctors

/*
 * WORK IN PROGRESS. USE WITH CAUTION
 *
 * TODO
 * - Currently limited to two inputs.
 * - Domain selection for secondary input?
 *
 * @Author Felix Maurer
 */

namespace coupling{
	template<unsigned int dim>
    class AsymmetricalFilterJunction;
}

template<unsigned int dim>
class coupling::AsymmetricalFilterJunction : public coupling::FilterSequence<dim> {
	public:
    	AsymmetricalFilterJunction( const coupling::IndexConversionMD2Macro<dim>* indexConversion,
						const tarch::utils::MultiMDService<dim>& multiMDService,
						const char* name,
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* > primaryInputCellVector, // primary input of sequence. 
						std::vector<tarch::la::Vector<dim, unsigned int>> primaryCellIndices, 
						const std::vector<coupling::datastructures::MacroscopicCell<dim>* >	secondaryInputCellVector, // additional data, presented as macro cells as well
						std::vector<tarch::la::Vector<dim, unsigned int>> secondaryCellIndices, 
						tarch::la::Vector<dim, unsigned int> primaryDomainStart,
						tarch::la::Vector<dim, unsigned int> primaryDomainEnd,
						std::array<bool, 7> filteredValues
		):
		coupling::FilterSequence<dim>(	indexConversion, 
										multiMDService, 
										name, 
										primaryInputCellVector,
									   	primaryCellIndices, 
										primaryDomainStart, 
										primaryDomainEnd, 
										filteredValues),
		_cellIndices_secondary(secondaryCellIndices),
		_inputCellVector_secondary(secondaryInputCellVector)
		{	
			#ifdef DEBUG_FILTER_JUNCTION_ASYM
        	std::cout << PRINT_PREFIX() << "Begin initialization." << std::endl;
        	#endif

			//allocate and init secondary cell vectors
			for(auto cell: _inputCellVector_secondary){
					_cellVector1_secondary.push_back(new coupling::datastructures::MacroscopicCell<dim>(*cell));
    				_cellVector2_secondary.push_back(new coupling::datastructures::MacroscopicCell<dim>(*cell));
			}
			#ifdef DEBUG_FILTER_JUNCTION_ASYM
				std::cout << PRINT_PREFIX() << "Initialized secondary cell vectors." << std::endl;
				std::cout << PRINT_PREFIX() << "First element of _cellVector1_secondary after init: " << _cellVector1_secondary[0] << std::endl;
				std::cout << PRINT_PREFIX() << "First element of _cellVector2_secondary after init: " << _cellVector2_secondary[0] << std::endl;
    		#endif 

			coupling::FilterSequence<dim>::_isModifiable = false; //Dynamic filters are not yet supported. TODO
    	}

    	~AsymmetricalFilterJunction(){
			for(auto secondarycell : _cellVector1_secondary) delete secondarycell;
			for(auto secondarycell : _cellVector2_secondary) delete secondarycell;
		}

		/*
		 * This function is very similar to the interface's. Check coupling::FilterSequence for more details.
		 */
		int loadFiltersFromXML(tinyxml2::XMLElement* sequenceNode) override;	

		void printFilters() override {
			std::cout << "Junctors in asymmetrical junction " << coupling::FilterSequence<dim>::_name << ": ";
			for(auto f : coupling::FilterSequence<dim>::_filters) std::cout << f->getType() << " ";
			std::cout << std::endl;
		}

		std::string PRINT_PREFIX() const override {
			return std::string("	AFJ(").std::string::append(coupling::FilterSequence<dim>::_name).std::string::append("): ");
		}


	private:
		//all secondary cell vectors use the same indexing
		std::vector<tarch::la::Vector<dim, unsigned int>> _cellIndices_secondary;

		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _inputCellVector_secondary;
		
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector1_secondary;
		std::vector<coupling::datastructures::MacroscopicCell<dim>* > _cellVector2_secondary;
};

//inlcude implementation
#include "coupling/filtering/sequencing/AsymmetricalFilterJunction.cpph"
