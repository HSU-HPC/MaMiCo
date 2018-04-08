// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_NIEVELOCITYIMPOSITION_H_
#define _MOLECULARDYNAMICS_COUPLING_NIEVELOCITYIMPOSITION_H_

#include "coupling/MomentumInsertion.h"
#include "coupling/cell-mappings/ComputeAvgForceAndVelocity.h"
#include "coupling/cell-mappings/NieVelocityImpositionMapping.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
  template<class LinkedCell,unsigned int dim>
  class NieVelocityImposition;
}


/** Velocity imposition scheme following the respective paper by Nie et al., J. Fluid. Mech. 500, 2004
 *
 *  @author Philipp Neumann
 */
template<class LinkedCell,unsigned int dim>
class coupling::NieVelocityImposition:
public coupling::MomentumInsertion<LinkedCell,dim> {
  public:
    NieVelocityImposition(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim>& indexConversion,
      const unsigned int &outermostLayer,
      const unsigned int &innermostLayer
    ):
    coupling::MomentumInsertion<LinkedCell,dim>(mdSolverInterface),
    _indexConversion(indexConversion),
    _outermostLayer(outermostLayer),
    _innermostLayer(innermostLayer){}
    virtual ~NieVelocityImposition(){}

    /** returns the number of MD steps between subsequent momentum insertions */
    virtual unsigned int getTimeIntervalPerMomentumInsertion() const { return 1;}

    virtual void insertMomentum(
      coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim>& cell,
      const unsigned int& currentLocalMacroscopicCellIndex
    ) const {
      // nop if this is not an imposition cell
      if (!isInsideImpositionLayer(currentLocalMacroscopicCellIndex)){ return;}
      // set continuum velocity
      tarch::la::Vector<dim,double> continuumVelocity(cell.getMicroscopicMomentum());

      coupling::cellmappings::ComputeAvgForceAndVelocity<LinkedCell,dim> computeForceAndVelocity(coupling::MomentumInsertion<LinkedCell,dim>::_mdSolverInterface);
      cell.iterateConstCells(computeForceAndVelocity);
      const tarch::la::Vector<dim,double> avgVel(computeForceAndVelocity.getAvgVelocity());
      const tarch::la::Vector<dim,double> avgF(computeForceAndVelocity.getAvgForce());

      coupling::cellmappings::NieVelocityImpositionMapping<LinkedCell,dim> velocityImposition(continuumVelocity,avgVel,avgF,coupling::MomentumInsertion<LinkedCell,dim>::_mdSolverInterface);
      cell.iterateCells(velocityImposition);
    }

  private:
    /** returns true if the local cell at index currentLocalMacroscopicCell is inside the layer of imposition cells, given by outermostLayer and innermostLayer. For, e.g.,
     *  outermostLayer=2 and innermostLayer=3, the layers for imposition are located in the 3rd and 4th strip of cells (we start counting from cell layer=0 which corresponds to the outermost,
     *  actually ghost-layer of cells which surrounds the MD domain).
     */
    bool isInsideImpositionLayer(const unsigned int &currentLocalMacroscopicCellIndex) const{
      const tarch::la::Vector<dim,unsigned int> globalNumberMacroscopicCells(_indexConversion.getGlobalNumberMacroscopicCells());
      const tarch::la::Vector<dim,unsigned int> globalCellIndex(_indexConversion.getGlobalVectorCellIndex(_indexConversion.convertLocalToGlobalCellIndex(currentLocalMacroscopicCellIndex)));
      bool isInsideLayer=false;
      // all directions
      for (unsigned int d = 0; d < dim; d++){
        isInsideLayer= isInsideLayer ||  ( (globalCellIndex[d]>= _outermostLayer) && (globalCellIndex[d]<=_innermostLayer) )
                      || ( (globalCellIndex[d]>= 1+globalNumberMacroscopicCells[d]-_innermostLayer) && (globalCellIndex[d] <= 1+globalNumberMacroscopicCells[d]-_outermostLayer) );
      }
      return isInsideLayer;
    }

    const coupling::IndexConversion<dim>& _indexConversion;
    const unsigned int _outermostLayer;
    const unsigned int _innermostLayer;
};


#endif // _MOLECULARDYNAMICS_COUPLING_VELOCITYGRADIENTRELAXATION_H_
