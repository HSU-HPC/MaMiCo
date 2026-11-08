// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
#define _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_

#include "coupling/interface/MoleculeIterator.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDMolecule.h"
#include "coupling/interface/impl/SimpleMD/SimpleMDLinkedCellWrapper.h"
#include "simplemd/services/MoleculeService.h"
#include "simplemd/LinkedCell.h"
#include "simplemd/Molecule.h"

namespace coupling {
namespace interface {

/** iterates over molecules in a SimpleMD linked cell.
 *  @author Philipp Neumann
 */
class SimpleMDMoleculeIterator : public MoleculeIterator<SimpleMDLinkedCellWrapper, MD_DIM> {
public:
  SimpleMDMoleculeIterator(SimpleMDLinkedCellWrapper& cellWrapper, simplemd::services::MoleculeService& moleculeService)
      : coupling::interface::MoleculeIterator<SimpleMDLinkedCellWrapper, MD_DIM>(cellWrapper), _moleculeService(moleculeService), _moleculeIndex(0),
        _buffer(NULL) {}
  virtual ~SimpleMDMoleculeIterator() {}

  /** sets the iterator to the first element */
  void begin() { _moleculeIndex = 0; }

  /** returns false, if the iterator reached the end of the molecule list */
  bool continueIteration() const {
    auto cell = _moleculeService.getContainer()[_cell.getCellIndex()];
    return (_moleculeIndex < cell.numMolecules());
  }

  /** sets the iterator to the next molecule */
  void next() {
    // check if we need to update molecule information and reset flag in this case
    _moleculeIndex++;
  }

  /** returns a reference to the molecule that this iterator currently points to */
  coupling::interface::Molecule<MD_DIM>& get() {
    auto cell = _moleculeService.getContainer()[_cell.getCellIndex()];
    auto molecule = &_moleculeService.getContainer().getMoleculeAt(cell.getIndex(), _moleculeIndex);
    _buffer.setMolecule(molecule);
    return _buffer;
  }

  const coupling::interface::Molecule<MD_DIM>& getConst() { return get(); }

private:
  simplemd::services::MoleculeService& _moleculeService;
  size_t _moleculeIndex;
  coupling::interface::SimpleMDMolecule _buffer;
};

} // namespace interface
} // namespace coupling

#endif // _MOLECULARDYNAMIC_COUPLING_INTERFACE_SIMPLEMDMOLECULEITERATOR_CPP_
