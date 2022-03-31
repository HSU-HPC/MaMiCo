// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_KINETICENERGYCONTROLLER_H_
#define _MOLECULARDYNAMICS_COUPLING_KINETICENERGYCONTROLLER_H_

#include "coupling/datastructures/MacroscopicCell.h"
#include "coupling/cell-mappings/ComputeKineticEnergyMapping.h"
#include "coupling/cell-mappings/SetKineticEnergyMapping.h"
#include "coupling/cell-mappings/ComputeMomentumMapping.h"
#include "coupling/cell-mappings/ComputeMassMapping.h"
#include "coupling/cell-mappings/ComputeTemperatureMapping.h"
#include "coupling/cell-mappings/SetTemperatureMapping.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class KineticEnergyController;
}

/** controles and regulates the kinetic energy of the MD system.
 *
 *  @author Philipp Neumann
 */
template <class LinkedCell, unsigned int dim>
class coupling::KineticEnergyController {
public:
  KineticEnergyController(coupling::interface::MDSolverInterface<
      LinkedCell, dim> *const mdSolverInterface)
      : _mdSolverInterface(mdSolverInterface) {}
  ~KineticEnergyController() {}

  /** computes and returns the kinetic energy within a macroscopic cell */
  double
  computeKineticEnergy(coupling::datastructures::MacroscopicCellWithLinkedCells<
      LinkedCell, dim> &cell) const {
    coupling::cellmappings::ComputeKineticEnergyMapping<LinkedCell, dim>
        computeKineticEnergyMapping(_mdSolverInterface);
    cell.iterateConstCells(computeKineticEnergyMapping);
    return computeKineticEnergyMapping.getKineticEnergy();
  }

  /** computes and returns the temperature within a macroscopic cell */
  double
  computeTemperature(coupling::datastructures::MacroscopicCellWithLinkedCells<
      LinkedCell, dim> &cell) const {
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
        computeMomentumMapping(_mdSolverInterface);
    cell.iterateConstCells(computeMomentumMapping);
    tarch::la::Vector<dim, double> meanVelocity =
        computeMomentumMapping.getMeanVelocity();

    coupling::cellmappings::ComputeTemperatureMapping<LinkedCell, dim>
        computeTemperatureMapping(meanVelocity, _mdSolverInterface);
    cell.iterateConstCells(computeTemperatureMapping);
    return computeTemperatureMapping.getTemperature();
  }

  /** sets the kinetic energy within a macroscopic cell to the value specified
   * by kineticEnergy.
   *  Therefore, the mean velocity is computed first. Afterwards the deviation
   * from the mean velocity
   *  are rescaled such that momentum is conserved.
   */
  void
  setKineticEnergy(coupling::datastructures::MacroscopicCellWithLinkedCells<
                       LinkedCell, dim> &cell,
                   const double &kineticEnergy) const {
    // determine mass, momentum and old kinetic energy
    coupling::cellmappings::ComputeMassMapping<LinkedCell, dim>
        computeMassMapping(_mdSolverInterface);
    cell.iterateConstCells(computeMassMapping);
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
        computeMomentumMapping(_mdSolverInterface);
    cell.iterateConstCells(computeMomentumMapping);
    coupling::cellmappings::ComputeKineticEnergyMapping<LinkedCell, dim>
        computeKineticEnergyMapping(_mdSolverInterface);
    cell.iterateConstCells(computeKineticEnergyMapping);

    // set new kinetic energy
    unsigned int numberParticles = computeMassMapping.getNumberOfParticles();
    tarch::la::Vector<dim, double> meanVelocity =
        computeMomentumMapping.getMeanVelocity();
    double oldKineticEnergy = computeKineticEnergyMapping.getKineticEnergy();
    coupling::cellmappings::SetKineticEnergyMapping<LinkedCell, dim>
        setKineticEnergyMapping(oldKineticEnergy, kineticEnergy,
                                numberParticles, meanVelocity);
    cell.iterateCells(setKineticEnergyMapping);
  }

  /** sets the temperature within the cell to 'temperature'. Here we just scale
   * the deviation from the mean velocity
   *  accordingly, i.e. we set: v_molecule = v_mean +
   * sqrt(temperature/current_temperature)*(v_molecule-v_mean)
   */
  void setTemperature(coupling::datastructures::MacroscopicCellWithLinkedCells<
                          LinkedCell, dim> &cell,
                      const double &temperature) const {
    coupling::cellmappings::ComputeMomentumMapping<LinkedCell, dim>
        computeMomentumMapping(_mdSolverInterface);
    cell.iterateConstCells(computeMomentumMapping);
    tarch::la::Vector<dim, double> meanVelocity =
        computeMomentumMapping.getMeanVelocity();

    coupling::cellmappings::ComputeTemperatureMapping<LinkedCell, dim>
        computeTemperatureMapping(meanVelocity, _mdSolverInterface);
    cell.iterateConstCells(computeTemperatureMapping);
    double currentTemperature = computeTemperatureMapping.getTemperature();

    coupling::cellmappings::SetTemperatureMapping<LinkedCell, dim>
        setTemperatureMapping(currentTemperature, temperature, meanVelocity,
                              _mdSolverInterface);
    cell.iterateCells(setTemperatureMapping);
  }

private:
  coupling::interface::MDSolverInterface<LinkedCell,
                                         dim> *const _mdSolverInterface;
};
#endif // _MOLECULARDYNAMICS_COUPLING_KINETICENERGYCONTROLLER_H_
