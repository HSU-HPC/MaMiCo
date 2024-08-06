// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_SERVICES_COUPLINGCELLSERVICEDUMMY_H_
#define _MOLECULARDYNAMICS_COUPLING_SERVICES_COUPLINGCELLSERVICEDUMMY_H_

#include "coupling/services/CouplingCellService.h"

#include "coupling/sendrecv/DataExchangeFromMD2Macro.h"
#include "coupling/sendrecv/DataExchangeFromMacro2MD.h"
#include "coupling/sendrecv/FromMD2MacroRecvOnly.h"
#include "coupling/sendrecv/FromMacro2MDSendOnly.h"

namespace coupling {
namespace services {
template <unsigned int dim> class CouplingCellServiceDummy;
}
} // namespace coupling

/** class for functionality of data exchange in hybrid Micro-Macro simulations.
 * This implementation assumes that MD does not run on this rank. We therefore
 * do not have any particle-related operations implemented here. Besides, no
 * grid for data exchange is allocated (except for some buffer cells inside the
 * data exchange objects from sendrecv).
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::services::CouplingCellServiceDummy : public coupling::services::CouplingCellService<dim> {
public:
  CouplingCellServiceDummy(unsigned int ID,
                           coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface, // interface to macroscopic solver
                           tarch::la::Vector<dim, unsigned int> numberProcesses,                             // number of processes in all directions
                           unsigned int rank,                                                                // current rank
                           tarch::la::Vector<dim, double> globalMDDomainSize,                                // domain size of MD simulation -> required for

                           tarch::la::Vector<dim, double> globalMDDomainOffset, // domain offset of MD simulation

                           const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration, // configuration for parallel topology
                           const coupling::configurations::CouplingCellConfiguration<dim>& couplingCellConfiguration,    // configuration for coupling cells
                                                                                                                         // and respective plotting
                           unsigned int topologyOffset)
      : coupling::services::CouplingCellService<dim>(ID), _macroscopicSolverInterface(macroscopicSolverInterface),
        _deFromMacro2MD(_macroscopicSolverInterface, ID), _deFromMD2Macro(_macroscopicSolverInterface, ID) {
    if (_macroscopicSolverInterface == NULL) {
      std::cout << "ERROR "
                   "coupling::services::CouplingCellServiceDummy::"
                   "CouplingCellServiceDummy(): "
                   "_macroscopicSolverInterface==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  virtual ~CouplingCellServiceDummy() {}

  void sendFromMacro2MD(const coupling::datastructures::FlexibleCellContainer<dim>& macro2MDBuffer) override {
    _fromMacro2MD.sendFromMacro2MD(_deFromMacro2MD, macro2MDBuffer);
  }
  double sendFromMD2Macro(const coupling::datastructures::FlexibleCellContainer<dim>& couplingCellContainerFromMacroscopicSolver) override {
    _fromMD2Macro.sendFromMD2Macro(_deFromMD2Macro, couplingCellContainerFromMacroscopicSolver);
    return 0;
  }
  double applyFilterPipeline() override { return 0; }
  void sendFromMacro2MDPreProcess() override {}
  void sendFromMacro2MDPostProcess() override {}
  void sendFromMD2MacroPreProcess() override {}
  void sendFromMD2MacroPostProcess() override {}
  void processInnerCouplingCellAfterMDTimestep() override {}
  void computeAndStoreTemperature(double temperature) override {}
  void applyTemperatureToMolecules(unsigned int t) override {}
  void distributeMass(unsigned int t) override {}
  void distributeMomentum(unsigned int t) override {}
  void perturbateVelocity() override {}
  void applyBoundaryForce(unsigned int t) override {}
  void plotEveryMicroscopicTimestep(unsigned int t) override {}
  void plotEveryMacroscopicTimestep(unsigned int t) override {}

private:
  /** interface for macroscopic solver */
  coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface;

  /** for quantity transfer between solvers */
  coupling::sendrecv::FromMacro2MDSendOnly<coupling::datastructures::CouplingCell<dim>, dim> _fromMacro2MD;
  coupling::sendrecv::DataExchangeFromMacro2MD<dim> _deFromMacro2MD;
  coupling::sendrecv::FromMD2MacroRecvOnly<coupling::datastructures::CouplingCell<dim>, dim> _fromMD2Macro;
  coupling::sendrecv::DataExchangeFromMD2Macro<dim> _deFromMD2Macro;
};

#endif // _MOLECULARDYNAMICS_COUPLING_COUPLINGCELLSERVICEDUMMY_H_
