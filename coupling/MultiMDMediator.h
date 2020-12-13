// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_MULTIMDMEDIATOR_H_
#define _COUPLING_MULTIMDMEDIATOR_H_

#include "coupling/InstanceHandling.h"
#include "coupling/services/MultiMDCellService.h"

namespace coupling {
  template<class LinkedCell, unsigned int> class MultiMDMediator;
}

/** Class to handle interaction between MultiMDCellService and InstanceHandling
 * This is currently mainly used for activating/deactivating and addtion/removal of md simulations
 */
template<class LinkedCell, unsigned int dim>
class coupling::MultiMDMediator {
public:

  MultiMDMediator(coupling::services::MultiMDCellService<LinkedCell, dim> & multiMDCellService, 
                  coupling::InstanceHandling<LinkedCell, dim> & instanceHandling, 
                  tarch::utils::MultiMDService<dim> & multiMDService,
                  coupling::interface::MacroscopicSolverInterface<dim> * macroscopicSolverInterface) 
    : _multiMDCellService(multiMDCellService), _instanceHandling(instanceHandling), _multiMDService(multiMDService),
      _listActiveMDSimulations(_multiMDService.getNumberLocalComms(), std::vector<bool>()),
      _macroscopicSolverInterface(macroscopicSolverInterface),
      _nextFreeBlock(_multiMDService.getNumberLocalComms()-1)
  {
    for(auto & group : _listActiveMDSimulations) {
      group = std::vector<bool>(_multiMDService.getLocalNumberOfMDSimulations(), true);
    }
  }


  ~MultiMDMediator() {
    if(_macroscopicSolverInterface != nullptr) {
      delete _macroscopicSolverInterface;
      _macroscopicSolverInterface = nullptr;
    }
  }


  /** Autmatically add another MDSimulation trying to keep the number of MD simulation across communicators balanced.
   */ 
  void addMDSimulation();


  /** Add MD simulation on specified communicator.
   * TODO How to handle already occupied locations??
   */
  void addMDSimulation(const unsigned int &);


  /** Add n MD simulations **/ //TODO
  void addNMDSimulations(const unsigned int &);


  /** ADD n MD Simulations to specified communicator */ //TODO
  void addNMDSimulations(const unsigned int &, const unsigned int &);


  /** Automatically remove MD simulation trying to keep the number of MD simulations across communicator balanced.
   */
  void rmMDSimulation();


  /** Remove MD simulation with specified global identifier.
   * TODO How to handle removal of empty places??
   */
  void rmMDSimulation(const unsigned int &, const unsigned int &);


  /** Remove MD Simulation on specific communicator
   */
  void rmMDSimulation(const unsigned int &);


  /** Remove N MD Simulations **/
  void rmNMDSimulations(const unsigned int &);

  /** Remove AL simulations on this communicator. 
   * This is the only way to do so!
  */
  void shutdownCommunicator(const unsigned int &);

private:


  /** Add one block of free simulations which is evenly sliced over communicator groups. */
  void addMDSimulationBlock();


  /** Remove one sliced block of simulations if all of them are freed. */
  void rmMDSimulationBlock();


  /** Reserve one slot in order to add a MD simulation.
   * This method reserves slots in a round robin fashion, starting at the communicator with the highest ID.
   */
  unsigned int reserveNextFreeSlot();


  /** Determine slot which has to be freed next in order to keep distribution of simulation as even as possible.
   */
  unsigned int getLastReservedSlot();


  /** Find global number of active simulations */
  unsigned int getNumberOfActiveMDSimulations();


  /** Find number of active simulations local to communicator */
  unsigned int getNumberOfActiveMDSimulations(const unsigned int);


  unsigned int getAvgNumberOfActiveMDSimulations();


  /** Try to find a communicator that has a relatively low number
   * of active MD simulations
   * This would be a communicator with a number of active simulations
   * lower than average. 
   * */
  unsigned int findFirstCommWithLowLoad();


  unsigned int findFirstCommWithHighLoad();


  /** On given communicator find an inactive index */
  unsigned int findInactiveLocalIndex(const unsigned int &);


  unsigned int findActiveLocalIndex(const unsigned int &);


  coupling::services::MultiMDCellService<LinkedCell, dim> & _multiMDCellService;
  coupling::InstanceHandling<LinkedCell, dim> & _instanceHandling;
  tarch::utils::MultiMDService<dim> & _multiMDService;
  std::vector<std::vector<bool> > _listActiveMDSimulations; // global list of active (true) and inactive (false) simulations
  coupling::interface::MacroscopicSolverInterface<dim> * _macroscopicSolverInterface;
  unsigned int _nextFreeBlock;
};

#include "coupling/MultiMDMediator.cpph"

#endif //_COUPLING_MULTIMDMEDIATOR_H_