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
      _macroscopicSolverInterface(macroscopicSolverInterface)
  {
    for(auto & group : _listActiveMDSimulations) {
      group = std::vector<bool>(_multiMDService.getLocalNumberOfMDSimulations(), true);
    }
  }

  /** Autmatically add another MDSimulation trying to keep the number of MD simulation across communicators balanced
   * by computing the average number of active simulations on communicators.
   */ 
  void addMDSimulation();


  /** Add MD simulation on specified communicator.
   */
  void addMDSimulation(const unsigned int &);


  /** Add n MD simulations **/
  void addNMDSimulations(const unsigned int &);


  /** ADD n MD Simulations to specified communicator */ //TODO
  void addNMDSimulations(const unsigned int &, const unsigned int &);


  /** Automatically remove MD simulation trying to keep the number of MD simulations across communicator balanced.
   */
  void rmMDSimulation();


  void forceRmMDSimulation(); 


  /** Remove MD simulation with specified global identifier.
   * TODO How to handle removal of empty places??
   */
  void rmMDSimulation(const unsigned int &, const unsigned int &);


  /** Remove MD Simulation on specific communicator
   */
  void rmMDSimulation(const unsigned int &);


  /** Remove N MD Simulations **/
  void rmNMDSimulations(const unsigned int &);


  void forceRmNMDSimulations(const unsigned int &);


  /** Remove ALL simulations on this communicator. 
   * This is the only way to do so!
  */
  void shutdownCommunicator(const unsigned int &);

private:


  /** Add one block of free simulations which is evenly sliced over communicator groups. */
  void addMDSimulationBlock();


  /** Remove one sliced block of simulations if all of them are freed. */
  void rmMDSimulationBlock();


  /** Find global number of active simulations */
  unsigned int getNumberOfActiveMDSimulations();


  /** Find number of active simulations local to communicator */
  unsigned int getNumberOfActiveMDSimulations(const unsigned int);


  /** Get the average number of active simulations running on the
   * communicator groups. The result will be correctly rounded to 
   * the nearest integer.
   */
  unsigned int getAvgNumberOfActiveMDSimulations();


  /** Try to find a communicator that has a relatively low number
   * of active MD simulations. This will be a communicator having
   * an average or below-average number of active simulations.
   * This method starts searching at the last communicator.
   * */
  unsigned int findFirstCommWithLowLoad();


  unsigned int findFirstCommWithHighLoad();


  /** On given communicator find an inactive index */
  unsigned int findInactiveLocalIndex(const unsigned int &);


  /** On a given communicator find a local index of an inactive md simulation */
  unsigned int findActiveLocalIndex(const unsigned int &);


  coupling::services::MultiMDCellService<LinkedCell, dim> & _multiMDCellService;
  coupling::InstanceHandling<LinkedCell, dim> & _instanceHandling;
  tarch::utils::MultiMDService<dim> & _multiMDService;
  std::vector<std::vector<bool> > _listActiveMDSimulations; // global list of active (true) and inactive (false) simulations
  coupling::interface::MacroscopicSolverInterface<dim> * _macroscopicSolverInterface;
};

#include "coupling/MultiMDMediator.cpph"

#endif //_COUPLING_MULTIMDMEDIATOR_H_