// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MoleculeIterator.h"
#include <iostream>

namespace coupling {
namespace cellmappings {
template <class LinkedCell, unsigned int dim> class MoleculeExtractor;
}
} // namespace coupling

/** extracts molecule information from a given macroscopic cell and stores all
 * molecule positions in a vector.
 *  This class is only meant for testing purposes!
 *  If you need individual access to molecules + some functionality, please do
 * so in a separate cell-mapping and
 *  call it from the respective transfer strategy instance.
 *	@brief This class extracts molecule information from a given macroscopic cell
 * and stores all molecule positions in a vector.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim> class coupling::cellmappings::MoleculeExtractor {
public:
  /** Constructor
   *	@param mdSolverInterface
   */
  MoleculeExtractor(coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface) : _mdSolverInterface(mdSolverInterface) {
    _molecules.clear();
  }

  /** Destructor */
  ~MoleculeExtractor() {}

  /** clear the vector of the molecule position, to make sure it is empty before
the iteration starts.
         */
  void beginCellIteration() { _molecules.clear(); }

  /** empty function
   */
  void endCellIteration() {}

  /** This function stores all molecule positions in a vector
   *	@param cell
   *	@param cellIndex
   */
  void handleCell(LinkedCell &cell, unsigned int &cellIndex) {
    coupling::interface::MoleculeIterator<LinkedCell, dim> *it = _mdSolverInterface->getMoleculeIterator(cell);
    it->begin();
    while (it->continueIteration()) {
      _molecules.push_back(it->getConst().getPosition());
      it->next();
    }
    delete it;
  }

  /** returns access to the extracted molecules.
   *	@return _molecules
   */
  const std::vector<tarch::la::Vector<dim, double>> &getExtractedMolecules() const { return _molecules; }

private:
  coupling::interface::MDSolverInterface<LinkedCell, dim> *const _mdSolverInterface;
  std::vector<tarch::la::Vector<dim, double>> _molecules;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_
