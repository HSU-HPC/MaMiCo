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
                                                                                                             // index conversion
                           tarch::la::Vector<dim, double> globalMDDomainOffset,                              // domain offset of MD simulation -> required
                                                                                                             // for index conversion
                           const coupling::configurations::ParallelTopologyConfiguration& parallelTopologyConfiguration, // configuration for parallel topology
                           const coupling::configurations::CouplingCellConfiguration<dim>& couplingCellConfiguration,    // configuration for coupling cells
                                                                                                                         // and respective plotting
                           unsigned int topologyOffset)
      : coupling::services::CouplingCellService<dim>(ID),
        // index conversion should be the very first thing getting initialised!
        _indexConversion(initIndexConversion(couplingCellConfiguration.getCouplingCellSize(), numberProcesses, rank, globalMDDomainSize, globalMDDomainOffset,
                                             parallelTopologyConfiguration.getParallelTopologyType(), topologyOffset)),
        _macroscopicSolverInterface(macroscopicSolverInterface), _deFromMacro2MD(_macroscopicSolverInterface, _indexConversion, ID),
        _deFromMD2Macro(_macroscopicSolverInterface, _indexConversion, ID) {
    if (_macroscopicSolverInterface == NULL) {
      std::cout << "ERROR "
                   "coupling::services::CouplingCellServiceDummy::"
                   "CouplingCellServiceDummy(): "
                   "_macroscopicSolverInterface==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  virtual ~CouplingCellServiceDummy() {
    if (_indexConversion != NULL) {
      delete _indexConversion;
      _indexConversion = NULL;
    }
  }

  virtual void sendFromMacro2MD(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver,
                                const unsigned int* const globalCellIndicesFromMacroscopicSolver) {
    _fromMacro2MD.sendFromMacro2MD(*_indexConversion, _deFromMacro2MD, couplingCellsFromMacroscopicSolver, globalCellIndicesFromMacroscopicSolver);
  }
  virtual double sendFromMD2Macro(const std::vector<coupling::datastructures::CouplingCell<dim>*>& couplingCellsFromMacroscopicSolver,
                                  const unsigned int* const globalCellIndicesFromMacroscopicSolver) {
    _fromMD2Macro.sendFromMD2Macro(*_indexConversion, _deFromMD2Macro, couplingCellsFromMacroscopicSolver, globalCellIndicesFromMacroscopicSolver);
    return 0;
  }
  virtual double applyFilterPipeline() { return 0; }
  virtual void sendFromMacro2MDPreProcess() {}
  virtual void sendFromMacro2MDPostProcess() {}
  virtual void sendFromMD2MacroPreProcess() {}
  virtual void sendFromMD2MacroPostProcess() {}
  virtual void processInnerCouplingCellAfterMDTimestep() {}
  virtual void computeAndStoreTemperature(double temperature) {}
  virtual void applyTemperatureToMolecules(unsigned int t) {}
  virtual void distributeMass(unsigned int t) {}
  virtual void distributeMomentum(unsigned int t) {}
  virtual void perturbateVelocity() {}
  virtual void applyBoundaryForce(unsigned int t) {}
  virtual void plotEveryMicroscopicTimestep(unsigned int t) {}
  virtual void plotEveryMacroscopicTimestep(unsigned int t) {}
  virtual const coupling::IndexConversion<dim>& getIndexConversion() const { return *_indexConversion; }

  virtual void updateIndexConversion(const unsigned int& topologyOffset) {
    auto* newIndexConversion = initIndexConversion(_indexConversion->getCouplingCellSize(), _indexConversion->getNumberProcesses(),
                                                   _indexConversion->getThisRank(), _indexConversion->getGlobalMDDomainSize(),
                                                   _indexConversion->getGlobalMDDomainOffset(), _indexConversion->getParallelTopologyType(), topologyOffset);

    delete _indexConversion;
    _indexConversion = newIndexConversion;

    _deFromMacro2MD.setIndexConversion(_indexConversion);
    _deFromMD2Macro.setIndexConversion(_indexConversion);
  }

private:
  /** this is currently a copy-paste version of CouplingCellService's
   * initIndexConversion */
  coupling::IndexConversion<dim>* initIndexConversion(tarch::la::Vector<dim, double> couplingCellSize, tarch::la::Vector<dim, unsigned int> numberProcesses,
                                                      unsigned int rank, tarch::la::Vector<dim, double> globalMDDomainSize,
                                                      tarch::la::Vector<dim, double> globalMDDomainOffset,
                                                      coupling::paralleltopology::ParallelTopologyType parallelTopologyType,
                                                      unsigned int topologyOffset) const {
    tarch::la::Vector<dim, unsigned int> globalNumberCouplingCells(0);
    for (unsigned int d = 0; d < dim; d++) {
      globalNumberCouplingCells[d] = (unsigned int)floor(globalMDDomainSize[d] / couplingCellSize[d] + 0.5);
      if (fabs(globalNumberCouplingCells[d] * couplingCellSize[d] - globalMDDomainSize[d]) > 1e-13) {
        std::cout << "coupling::services::CouplingCellServiceDummy::"
                     "initIndexConversion(): Deviation of domain size > 1e-13!"
                  << std::endl;
      }
    }
    coupling::IndexConversion<dim>* ic = new coupling::IndexConversion<dim>(globalNumberCouplingCells, numberProcesses, rank, globalMDDomainSize,
                                                                            globalMDDomainOffset, parallelTopologyType, topologyOffset);
    if (ic == NULL) {
      std::cout << "coupling::services::CouplingCellServiceImpl::"
                   "initIndexConversion(): ic==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return ic;
  }

  /** for index conversion */
  coupling::IndexConversion<dim>* _indexConversion;

  /** interface for macroscopic solver */
  coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface;

  /** for quantity transfer between solvers */
  coupling::sendrecv::FromMacro2MDSendOnly<coupling::datastructures::CouplingCell<dim>, dim> _fromMacro2MD;
  coupling::sendrecv::DataExchangeFromMacro2MD<dim> _deFromMacro2MD;
  coupling::sendrecv::FromMD2MacroRecvOnly<coupling::datastructures::CouplingCell<dim>, dim> _fromMD2Macro;
  coupling::sendrecv::DataExchangeFromMD2Macro<dim> _deFromMD2Macro;
};

#endif // _MOLECULARDYNAMICS_COUPLING_COUPLINGCELLSERVICEDUMMY_H_
