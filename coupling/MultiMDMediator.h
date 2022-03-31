// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUPLING_MULTIMDMEDIATOR_H_
#define _COUPLING_MULTIMDMEDIATOR_H_

#include "coupling/InstanceHandling.h"
#include "coupling/services/MultiMDCellService.h"

namespace coupling {
template <class LinkedCell, unsigned int> class MultiMDMediator;
}

/**
 *	@brief Class to handle interaction between MultiMDCellService and
 * InstanceHandling
 * 	This is currently mainly used for activating/deactivating and
 * addtion/removal of md simulations. Works and interacts with the class
 * coupling::InstanceHandling closely.
 * 	@tparam LinkedCell type of the cell
 * 	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *	@sa see also coupling::InstanceHandling
 *  @author Niklas Wittmer
 */
template <class LinkedCell, unsigned int dim> class coupling::MultiMDMediator {
public:
  MultiMDMediator(coupling::services::MultiMDCellService<LinkedCell, dim> &multiMDCellService, coupling::InstanceHandling<LinkedCell, dim> &instanceHandling,
                  tarch::utils::MultiMDService<dim> &multiMDService, coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface)
      : _multiMDCellService(multiMDCellService), _instanceHandling(instanceHandling), _multiMDService(multiMDService),
        _listActiveMDSimulations(_multiMDService.getNumberLocalComms(), std::vector<bool>()), _macroscopicSolverInterface(macroscopicSolverInterface) {
    for (auto &group : _listActiveMDSimulations) {
      group = std::vector<bool>(_multiMDService.getLocalNumberOfMDSimulations(), true);
    }
  }

  /** Autmatically add one MDSimulation trying to keep the number of MD
   * simulation across communicators balanced
   *  by computing the average number of active simulations on communicators.
   */
  void addMDSimulation();

  /**	Add one MD simulation on specified communicator.
   *	@param communicator
   */
  void addMDSimulation(const unsigned int &communicator);

  /** Add n MD simulations
   *  @param N number of new MD
   */
  void addNMDSimulations(const unsigned int &N);

  /** ADD n MD Simulations to specified communicator
   *  @param communicator
   *  @param N number of new MD
   */
  void addNMDSimulations(const unsigned int &communicator, const unsigned int &N);

  /** Automatically remove one MD simulation trying to keep the number of MD
   * simulations across communicator balanced.
   */
  void rmMDSimulation();

  /** Remove MD simulation with specified global identifier.
   *  @param communicator
   *  @param index global identifier
   *  @todo How to handle removal of empty places?
   */
  void rmMDSimulation(const unsigned int &communicator, const unsigned int &index);

  /** Remove MD Simulation on specific communicator
   *  @param communicator
   */
  void rmMDSimulation(const unsigned int &communicator);

  /** Remove N MD Simulations
   *  @param N number of new MD
   */
  void rmNMDSimulations(const unsigned int &N);

  /** Remove ALL simulations on a communicator.
   *  @param communicator
   */
  void shutdownCommunicator(const unsigned int &communicator);

  /** Find number of active simulations local to communicator
   *  @param communicator
   *  @return number of active simulations local to communicator
   */
  unsigned int getNumberOfActiveMDSimulations(const unsigned int communicator);

  /** Find global number of active simulations
   *  @return total number of active simulation
   */
  unsigned int getNumberOfActiveMDSimulations();

private:
  /** Add one block of free simulations which is evenly sliced over communicator
   * groups. */
  void addMDSimulationBlock();

  /** Remove one sliced block of simulations if all of them are freed. */
  void rmMDSimulationBlock();

  /** Get the average number of active simulations running on the
   * communicator groups. The result will be correctly rounded to
   * the nearest integer.
   *  @return average number of active simulations
   */
  unsigned int getAvgNumberOfActiveMDSimulations();

  /** Try to find a communicator that has a relatively low number
   * of active MD simulations. This will be a communicator having
   * an average or below-average number of active simulations.
   * This method starts searching at the last communicator.
   *  @return communicator with a relatively low number of active MD simulations
   */
  unsigned int findFirstCommWithLowLoad();

  /** Try to find a communicator that has a relatively high number
   * of active MD simulations. The average ov this communicator should
   * should have an average or higher-than-average number of active simulations.
   * This method begins lookup at the first communicator.
   *  @return communicator with a relatively high number of active MD
   * simulations
   */
  unsigned int findFirstCommWithHighLoad();

  /** On given communicator find an inactive index
   *  @param communicator
   *  @return inactive index
   */
  unsigned int findInactiveLocalIndex(const unsigned int &communicator);

  /** On a given communicator find a local index of an inactive md simulation
   *  @param communicator
   *  @return local index of an inactive md simulation
   */
  unsigned int findActiveLocalIndex(const unsigned int &communicator);

  coupling::services::MultiMDCellService<LinkedCell, dim> &_multiMDCellService;
  coupling::InstanceHandling<LinkedCell, dim> &_instanceHandling;
  tarch::utils::MultiMDService<dim> &_multiMDService;
  std::vector<std::vector<bool>> _listActiveMDSimulations; // global list of active (true) and inactive
                                                           // (false) simulations
  coupling::interface::MacroscopicSolverInterface<dim> *_macroscopicSolverInterface;
};

#include "coupling/MultiMDMediator.cpph"

#endif //_COUPLING_MULTIMDMEDIATOR_H_