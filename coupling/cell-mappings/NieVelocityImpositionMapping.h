// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_NIEVELOCITYIMPOSITIONMAPPING_H_
#define _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_NIEVELOCITYIMPOSITIONMAPPING_H_

#include <iostream>
#include "tarch/la/Matrix.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/IndexConversion.h"
#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/interface/Molecule.h"
#include "coupling/CouplingMDDefinitions.h"

namespace coupling {
  namespace cellmappings {
    template<class LinkedCell,unsigned int dim>
    class NieVelocityImpositionMapping;
  }
}


/** employs an acceleration based on velocity gradients (in time) using the forcing term of the molecules.
 *  We currently expect that valid force entries are provided in each molecule, that is all molecule-molecule
 *  interactions have previously been computed.
 *  The scheme follows the descriptions in the paper by Nie et al., J. Fluid. Mech. 500, 55-64, 2004.
 *  However, we only average over single grid cells (no averaging over one (periodic) dimension or similar tricks ;-) ).
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::cellmappings::NieVelocityImpositionMapping {
  public:
    NieVelocityImpositionMapping(
      const tarch::la::Vector<dim,double> &continuumVelocity,
      const tarch::la::Vector<dim,double> &avgMDVelocity,
      const tarch::la::Vector<dim,double> &avgForce,
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface
    ): _mdSolverInterface(mdSolverInterface),
       _continuumVelocity(continuumVelocity),
       _avgMDVelocity(avgMDVelocity),
       _avgForce(avgForce)
    {}


    virtual ~NieVelocityImpositionMapping(){}

    void beginCellIteration(){}

    void endCellIteration(){}

    void handleCell(LinkedCell& cell,const unsigned int &cellIndex){
      coupling::interface::MoleculeIterator<LinkedCell,dim> *it = _mdSolverInterface->getMoleculeIterator(cell);

      // compute average force contribution inside this linked cell
      it->begin();
      while(it->continueIteration()){
        coupling::interface::Molecule<dim> &wrapper(it->get());
        tarch::la::Vector<dim,double> force (wrapper.getForce());
        force = force - _avgForce - (_mdSolverInterface->getMoleculeMass()/_mdSolverInterface->getDt()) * (_avgMDVelocity - _continuumVelocity);
        wrapper.setForce(force);
        it->next();
      }
      delete it;
    }

  private:
    // MD solver interface, required for molecule iterator and molecule mass
    coupling::interface::MDSolverInterface<LinkedCell,dim> * const _mdSolverInterface;
    // current velocity in this macroscopic cell (=velocity from continuum solver)
    const tarch::la::Vector<dim,double> _continuumVelocity;
    // current avg. velocity of molecules
    const tarch::la::Vector<dim,double> _avgMDVelocity;
    // average force within this macroscopic cell
    const tarch::la::Vector<dim,double> _avgForce;
};


#endif // _MOLECULARDYNAMICS_COUPLING_CELLMAPPINGS_NIEVELOCITYIMPOSITIONMAPPING_H_
