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
	}
}

/*
 * TODO: is it really that smart to design this as a class and not just as a function?
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
						const coupling::interface::MacroscopicSolverInterface<dim> *msi);

	private:
		const simplemd::configurations::MolecularDynamicsConfiguration _simpleMDConfig;
		const coupling::configurations::MaMiCoConfiguration<dim> _mamicoConfig;
		const coupling::interface::MacroscopicSolverInterface<dim> *_msi; 
};
