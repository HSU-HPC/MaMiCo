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
template <class LinkedCell, unsigned int dim> class NieVelocityImposition;
}

/** @brief Velocity imposition scheme following the respective paper by Nie et
 * al., J. Fluid. Mech. 500, 2004
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 */
template <class LinkedCell, unsigned int dim> class coupling::NieVelocityImposition : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param mdSolverInterface interface to the md solver
   *  @param outermostLayer the index of the outermost cell layer
   *  @param innermostLayer the index of the innermost cell layer */
  NieVelocityImposition(coupling::interface::MDSolverInterface<LinkedCell, dim>* const mdSolverInterface, const unsigned int& outermostLayer,
                        const unsigned int& innermostLayer)
      : coupling::MomentumInsertion<LinkedCell, dim>(mdSolverInterface), _outermostLayer(outermostLayer), _innermostLayer(innermostLayer),
        _enableInnerImposition(false) {}

  /** @brief a simple destructor */
  virtual ~NieVelocityImposition() {}

  /** @brief momentum shall be inserted in every md time step, so this returns 1
   *  @returns the time step interval for momentum insertion, always 1 */
  unsigned int getTimeIntervalPerMomentumInsertion() const override { return 1; }

  /** @brief inserts momentum to a cell
   *  @param cell to the coupling cell will the momentum be inserted
   *  @param idx local linearised index for the
   * coupling cell */
  void insertMomentum(coupling::datastructures::CouplingCellWithLinkedCells<LinkedCell, dim>& cell, I02 idx) const override {
    // nop if this is not an imposition cell
    if (!isInsideImpositionLayer(idx)) {
      return;
    }
    // set continuum velocity
    tarch::la::Vector<dim, double> continuumVelocity(cell.getMicroscopicMomentum());

    coupling::cellmappings::ComputeAvgForceAndVelocity<LinkedCell, dim> computeForceAndVelocity(
        coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    cell.iterateConstCells(computeForceAndVelocity);
    const tarch::la::Vector<dim, double> avgVel(computeForceAndVelocity.getAvgVelocity());
    const tarch::la::Vector<dim, double> avgF(computeForceAndVelocity.getAvgForce());

    coupling::cellmappings::NieVelocityImpositionMapping<LinkedCell, dim> velocityImposition(continuumVelocity, avgVel, avgF,
                                                                                             coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    cell.iterateCells(velocityImposition);
  }

  void setInnerImposition(bool enable) override { _enableInnerImposition = enable; }

private:
  /** returns true if the local cell at index currentLocalCouplingCellIndex is
   * inside the layer of imposition cells, given by outermostLayer and
   * innermostLayer. For, e.g., outermostLayer=2 and innermostLayer=3, the
   * layers for imposition are located in the 3rd and 4th strip of cells (we
   * start counting from cell layer=0 which corresponds to the outermost,
   * actually ghost-layer of cells which surrounds the MD domain).
   *  @brief based on the cell index, the function tells if the cell is inside
   * the imposition layer
   *  @param globalCellIndex global linearised index of a coupling
   * cell to check
   *  @returns a bool, that indicates if the given cell index is located in the
   * imposition layer (true) or not (false) */
  bool isInsideImpositionLayer(I01 globalCellIndex) const {
    bool inner = true;
    tarch::la::Vector<dim, unsigned int> globalIndexUnsigned{globalCellIndex.get()};
    for (unsigned int d = 0; d < dim; d++)
      inner = inner && (globalIndexUnsigned[d] > _innermostLayer && globalIndexUnsigned[d] < 1 + I09::numberCellsInDomain[d] - _innermostLayer);
    bool outer = false;
    for (unsigned int d = 0; d < dim; d++)
      outer = outer || (globalIndexUnsigned[d] < _outermostLayer || globalIndexUnsigned[d] > 1 + I09::numberCellsInDomain[d] - _outermostLayer);
    return (_enableInnerImposition || !inner) && !outer;
  }

  /** @brief the index of the outermost cell layer*/
  const unsigned int _outermostLayer;
  /** @brief the index of the innermost cell layer*/
  const unsigned int _innermostLayer;
  bool _enableInnerImposition;
};

#endif // _MOLECULARDYNAMICS_COUPLING_VELOCITYGRADIENTRELAXATION_H_
