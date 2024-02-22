#pragma once

#include "coupling/services/CouplingCellService.h"

namespace coupling {
namespace interface {
/**
 *	@brief generic interface class for different microscopic (MD) solvers.
 *  @author Philipp Neumann
 */
class MDSimulation {
public:
  /** Destructor */
  virtual ~MDSimulation() {}

  /** switches coupling off*/
  virtual void switchOffCoupling() = 0;

  /** switches coupling on*/
  virtual void switchOnCoupling() = 0;

  /** simulates numberTimesteps time steps and starts at time step no.
   *firstTimestep
   *	@param numberTimesteps
   *	@param firstTimestep
   */
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) = 0;

  // simulates a single time step
  // virtual void simulateTimestep(const unsigned int &thisTimestep ){const
  // unsigned int steps=1; simulateTimesteps(thisTimestep,steps);} TODO BUG

  /** sortMoleculesIntoCells*/
  virtual void sortMoleculesIntoCells() = 0;

  /** setCouplingCellService
   *	@param couplingCellService
   */
  virtual void setCouplingCellService(coupling::services::CouplingCellService<MDSIMULATIONFACTORY_DIMENSION>* couplingCellService) = 0;

  /** initialises the _molecularDynamicsSimulation solver
   *	@sa simplemd::MolecularDynamicsSimulation::initServices()
   *	@todo Philipp ??
   */
  virtual void init() = 0;

  /** initialises the _molecularDynamicsSimulation solver
   *	@param multiMDService
   *	@param localMDSimulation
   *	@sa simplemd::MolecularDynamicsSimulation::initServices(const
   *tarch::utils::MultiMDService<MD_DIM>& multiMDService,unsigned int
   *localMDSimulation)
   *	@todo Philipp ??
   */
  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) = 0;

  /** shuts down the MD simulation*/
  virtual void shutdown() = 0;

  /** Saves the simulation result as check point in the file filestem
   *	@param filestem
   *	@param t
   */
  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) = 0;
};
} // namespace interface
} // namespace coupling
