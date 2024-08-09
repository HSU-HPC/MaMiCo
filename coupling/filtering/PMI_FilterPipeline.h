#pragma once

#include "coupling/filtering/FilterPipeline.h"

namespace coupling {
namespace filtering {
template <unsigned int dim> class PMI_FilterPipeline;
}
}

// Post Multi Instance FilterPipeline
template <unsigned int dim> class coupling::filtering::PMI_FilterPipeline : public coupling::filtering::FilterPipeline<coupling::datastructures::FlexibleCellContainer<dim>,dim> {

	PMI_FilterPipeline(const coupling::datastructures::FlexibleCellContainer<dim> inputCells, const tarch::utils::MultiMDService<dim>& multiMDService,
        const char* cfgpath) : FilterPipeline<coupling::datastructures::FlexibleCellContainer<dim>,dim>(multiMDService, cfgpath) {
		this->_md2MacroCells = inputCells;
		loadSequencesFromXML(this->_config.root->FirstChildElement("filter-pipeline")->FirstChildElement("post-multi-instance"));
	}
};
