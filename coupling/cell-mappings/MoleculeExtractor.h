// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_
#include <iostream>
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MoleculeIterator.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class MoleculeExtractor;
  }
}



/** extracts molecule information from a given macroscopic cell and stores all molecule positions in a vector.
 *  This class is only meant for testing purposes!
 *  If you need individual access to molecules + some functionality, please do so in a separate cell-mapping and
 *  call it from the respective transfer strategy instance.
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::MoleculeExtractor {
  public:
    MoleculeExtractor(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface):
    _mdSolverInterface(mdSolverInterface){ _molecules.clear(); }
    ~MoleculeExtractor(){}

    void beginCellIteration(){ _molecules.clear(); }
    void endCellIteration(){}
    void handleCell(LinkedCell &cell, unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim>* it = _mdSolverInterface->getMoleculeIterator(cell);
      it->begin();
      while(it->continueIteration()){
        _molecules.push_back(it->getConst().getPosition());
        it->next();
      }
      delete it;
    }

    /** returns access to the extracted molecules. */
    const std::vector<tarch::la::Vector<dim,double> >& getExtractedMolecules() const{return _molecules; }

  private:
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    std::vector<tarch::la::Vector<dim,double> > _molecules;
};
 

#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_MOLECULEEXTRACTOR_H_

