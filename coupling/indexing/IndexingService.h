#pragma once

#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MacroscopicSolverInterface.h"

//Include CellIndex template class definition
#include "CellIndex.h"

//Include non-member functions operating on indexes
#include "Operations.h"

namespace coupling {
	namespace indexing {

		template<unsigned int dim>
		class IndexingService;

		template<unsigned int dim>
		std::vector<unsigned int> 
		getRanksForGlobalIndex(const CellIndex<dim, BaseIndexType> &globalCellIndex, const tarch::la::Vector<dim, unsigned int> &globalNumberMacroscopicCells);

	}
}

/*
 * TODO: redesign as function
 *
 * TODO: check includes
 *
 * TODO: comment
 *
 * @author Felix Maurer
 */
template<unsigned int dim>
class coupling::indexing::IndexingService{
	public:
		IndexingService(const simplemd::configurations::MolecularDynamicsConfiguration &simpleMDConfig,
						const coupling::configurations::MaMiCoConfiguration<dim> &mamicoConfig,
						coupling::interface::MacroscopicSolverInterface<dim> *msi,
						const unsigned int rank);

	private:
		const simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
		const coupling::configurations::MaMiCoConfiguration<dim> _mamicoConfig;
		coupling::interface::MacroscopicSolverInterface<dim> *_msi; 
		const unsigned int _rank;
};
