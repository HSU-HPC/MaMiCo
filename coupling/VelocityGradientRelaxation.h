// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_VELOCITYGRADIENTRELAXATION_H_
#define _MOLECULARDYNAMICS_COUPLING_VELOCITYGRADIENTRELAXATION_H_

#include "coupling/MomentumInsertion.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/VelocityGradientRelaxationMapping.h"
#include "coupling/interface/MDSolverInterface.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class VelocityGradientRelaxation;
template <class LinkedCell, unsigned int dim> class VelocityGradientRelaxationTopOnly;
} // namespace coupling

/** In this particular case, however, the velocity is relaxed within a one
 * cell-wide strip around the molecular domain. For this purpose, the velocity
 * that shall be imposed in average in the cells that are within a three
 * cell-wide strip need to be known and stored in the
 *  microscopicMomentum-buffers (2 cells overlap with the MD simulation, and 1
 * more (ghost) cell layer surrounding the MD domain). The procedure then only
 * considers molecules that are located between the midpoints of the cells which
 * are in the two cell-wide boundary strip. For all these molecules, a
 * second-order interpolation of the average velocity at their respective
 * position is carried out and the molecules are then relaxed towards this
 * particular velocity.
 *  @brief carries out velocity relaxation (similar to the SetMomentumMapping
 * procedure).
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3 */
template <class LinkedCell, unsigned int dim> class coupling::VelocityGradientRelaxation : public coupling::MomentumInsertion<LinkedCell, dim> {
public:
  /** @brief a simple constructor
   *  @param relaxationParam defines the strength of the relaxation; 1 means the
   * molecules velocity is set to the new velocity; O.5 -> the velocity will be
   * set to the avaerage of the old and the new velocity
   *  @param mdSolverInterface interface for the md solver
   *  @param indexConversion instance of the indexConversion for the current
   * simulation
   *  @param macroscopicCells the macroscopic cells to apply the velocity
   * gradient relaxation */
  VelocityGradientRelaxation(double relaxationParam, coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface,
                             const coupling::IndexConversion<dim> &indexConversion,
                             const coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> *const macroscopicCells)
      : coupling::MomentumInsertion<LinkedCell, dim>(mdSolverInterface), _indexConversion(indexConversion), _macroscopicCells(macroscopicCells),
        _relaxationParam(relaxationParam) {}
  /** @brief a dummy destructor */
  virtual ~VelocityGradientRelaxation() {}

  /** @brief returns the time step interval for the momentum insertion, always
   * one for this method
   *  @return the time step interval for momentum insertion */
  virtual unsigned int getTimeIntervalPerMomentumInsertion() const { return 1; }

  /** This method does not conserve the kinetic energy of the respective
   * macroscopic cell. To conserve the energy as well, see the description of
   * MomentumController on details how to do that.
   *  @brief insertes the momentum to the cells according to the params and
   * velocity gradient relaxation method
   *  @param cell the macroscopic cell to insert momentum to
   *  @param currentLocalCellIndex the macroscopic cell's index */
  virtual void insertMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                              const unsigned int &currentLocalMacroscopicCellIndex) const {
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> momentumMapping(coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    tarch::la::Vector<dim, double> oldVelocity(0.0);
    cell.iterateConstCells(momentumMapping);
    // set current average velocity within this cell
    oldVelocity = momentumMapping.getMeanVelocity();
    // set new momentum (based on velocity stored in microscopic
    // momentum-buffer)
    coupling::cellmappings::VelocityGradientRelaxationMapping<LinkedCell, dim> velocityGradientRelaxation(
        _relaxationParam, oldVelocity, currentLocalMacroscopicCellIndex, coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface, _indexConversion,
        _macroscopicCells);
    cell.iterateCells(velocityGradientRelaxation);
  }

protected:
  /** instance of the indexConversion for the current simulation */
  const coupling::IndexConversion<dim> &_indexConversion;
  /** the macroscopic cells to apply the velocity gradient relaxation */
  const coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> *const _macroscopicCells;
  /** defines the strength of the relaxation; 1 means the molecules velocity is
   * set to the new velocity;
   *                         O.5 -> the velocity will be set to the avaerage of
   * the old and the new velocity*/
  const double _relaxationParam;
};

/** In this particular case, however, the velocity is relaxed within a one
 * cell-wide strip around the molecular domain. For this purpose, the velocity
 * that shall be imposed in average in the cells that are within a three
 * cell-wide strip need to be known and stored in the
 *  microscopicMomentum-buffers (2 cells overlap with the MD simulation, and 1
 * more (ghost) cell layer surrounding the MD domain). The procedure then only
 * considers molecules that are located between the midpoints of the cells which
 * are in the two cell-wide boundary strip. For all these molecules, a
 * second-order interpolation of the average velocity at their respective
 * position is carried out and the molecules are then relaxed towards this
 * particular velocity.
 *  @brief carries out velocity relaxation (similar to the SetMomentumMapping
 * procedure).
 *  @author Philipp Neumann
 *  @tparam LinkedCell the LinkedCell class is given by the implementation of
 * linked cells in the molecular dynamics simulation
 *  @tparam dim  refers to the spacial dimension of the simulation, can be 1, 2,
 * or 3
 *  @todo Don't know what this does, just copied the comments from above*/
template <class LinkedCell, unsigned int dim> class coupling::VelocityGradientRelaxationTopOnly : public coupling::VelocityGradientRelaxation<LinkedCell, dim> {
public:
  /** @brief a simple constructor*/
  VelocityGradientRelaxationTopOnly(double relaxationParam, coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface,
                                    const coupling::IndexConversion<dim> &indexConversion,
                                    const coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> *const macroscopicCells)
      : coupling::VelocityGradientRelaxation<LinkedCell, dim>(relaxationParam, mdSolverInterface, indexConversion, macroscopicCells) {}

  /** @brief a dummy destructor */
  virtual ~VelocityGradientRelaxationTopOnly() {}

  /** This method does not conserve the kinetic energy of the respective
   * macroscopic cell. To conserve the energy as well, see the description of
   * MomentumController on details how to do that.
   *  @brief insertes the momentum to the cells according to the params and
   * velocity gradient relaxation method
   *  @param cell the macroscopic cell to insert momentum to
   *  @param currentLocalCellIndex the macroscopic cell's index */
  virtual void insertMomentum(coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell, dim> &cell,
                              const unsigned int &currentLocalMacroscopicCellIndex) const {
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim> momentumMapping(coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface);
    tarch::la::Vector<dim, double> oldVelocity(0.0);
    cell.iterateConstCells(momentumMapping);
    // set current average velocity within this cell
    oldVelocity = momentumMapping.getMeanVelocity();

    // set new momentum (based on velocity stored in microscopic
    // momentum-buffer)
    coupling::cellmappings::VelocityGradientRelaxationTopOnlyMapping<LinkedCell, dim> velocityGradientRelaxation(
        coupling::VelocityGradientRelaxation<LinkedCell, dim>::_relaxationParam, oldVelocity, currentLocalMacroscopicCellIndex,
        coupling::MomentumInsertion<LinkedCell, dim>::_mdSolverInterface, coupling::VelocityGradientRelaxation<LinkedCell, dim>::_indexConversion,
        coupling::VelocityGradientRelaxation<LinkedCell, dim>::_macroscopicCells);
    cell.iterateCells(velocityGradientRelaxation);
  }
};
#endif // _MOLECULARDYNAMICS_COUPLING_VELOCITYGRADIENTRELAXATION_H_
