#pragma once

#include "coupling/filtering/FilterPipeline.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class PI_FilterPipeline;
}
}

// Per Instance FilterPipeline
template <unsigned int dim> class coupling::filtering::PI_FilterPipeline : public coupling::filtering::FilterPipeline<coupling::datastructures::CellContainer<I14, dim>,dim> {

	PI_FilterPipeline(const coupling::datastructures::CellContainer<I02,dim> inputCells, const tarch::utils::MultiMDService<dim>& multiMDService,
        const char* cfgpath) : FilterPipeline<coupling::datastructures::CellContainer<I14, dim>,dim>(multiMDService, cfgpath) {

		I01 idx;
		coupling::datastructures::CouplingCell<3>* cell;
		for (auto pair : inputCells) {
		  std::tie(cell, idx) = pair;
		  if (I14::contains(idx))
		    this->_md2MacroCells << cell;
		  else
		    this->_outerCells << pair;
		}

		loadSequencesFromXML(this->_config.root->FirstChildElement("filter-pipeline")->FirstChildElement("per-instance"));
	}
};
