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

template<class LinkedCell, unsigned int dim>
class coupling::MultiMDMediator {
public:

  MultiMDMediator(coupling::services::MultiMDCellService<LinkedCell, dim> & multiMDCellService, 
                  coupling::InstanceHandling<dim> & instanceHandling, 
                  tarch::utils::MultiMDService<dim> & multiMDService) 
    : _multiMDCellService(multiMDCellService), _instanceHandling(instanceHandling), _multiMDService(multiMDService),
      _listActiveMDSimulations(_multiMDService.getNumberLocalComms(), std::vector<bool>()),
      _nextFreeBlock(_multiMDService.getNumberLocalComms())
    {
      for(auto & group : _listActiveMDSimulations) {
        group = std::vector<bool>(_multiMDService.getNumberLocalComms(), true);
      }
    }


  /** Autmatically add another MDSimulation trying to keep the number of MD simulation across communicators balanced.
   */ 
  void addMDSimulation(coupling::interface::MacroscopicSolverInterface<dim> *);


  /** Add MD simulation on specified communicator.
   * TODO How to handle already occupied locations??
   */
  void addMDSimulation(coupling::interface::MacroscopicSolverInterface<dim> *, const unsigned int &);


  /** Automatically remove MD simulation trying to keep the number of MD simulations across communicator balanced.
   */
  void rmMDSimulation();


  /** Remove MD simulation with specified global identifier.
   * TODO How to handle removal of empty places??
   */
  void rmMDSimulation(const unsigned int &, const unsigned int &);


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


private:
  coupling::services::MultiMDCellService<LinkedCell, dim> & _multiMDCellService;
  coupling::InstanceHandling<dim> & _instanceHandling;
  tarch::utils::MultiMDService<dim> & _multiMDService;
  std::vector<std::vector<bool> > _listActiveMDSimulations; // global list of active (true) and inactive (false) simulations
  unsigned int _nextFreeBlock;
};

#include "coupling/MultiMDMediator.cpph"

#endif //_COUPLING_MULTIMDMEDIATOR_H_