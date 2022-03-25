// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _INSTANCE_HANDLING_H_
#define _INSTANCE_HANDLING_H_

#include "coupling/interface/MDSimulationFactory.h"
#include "coupling/services/MultiMDCellService.h"
#include "tarch/utils/MultiMDService.h"

namespace coupling {
template <class LinkedCell, unsigned int dim> class InstanceHandling;
}

/** holds one vector of MDSimulation and one vector for MDSolverInterface.
 *Initialization, execution of MD time steps and shutdown are abstracted into
 *this class. In order to launch a new MF simulatio, a slot has to be chosen
 *first (either manualy or using coupling::MultiMDMediator). Then
 *MultiMDCellService initializes a new MacroscopicCellService via
 *MultiMDMediator. In the activated slot, a new MD simulation is launched. The
 *new MD instance has to be equilibrated first and then it can be coupled to the
 *simulation. In order to remove a MD simulation, the MD simulation and its
 *corresponding MDSolverInterface are shut down. Then, the respective instance
 *of the MacroscopicCellService is removed. Finally, the selected slot will be
 *set to inactive. This slot is now available again for the launch of a new MD
 *instance in the future.
 *	@brief Simulation slots are managed (i.e., added/removed) via this
 *class. Works and interacts with the class coupling::MultiMDMediator closely.
 * 	@tparam LinkedCell type of the cell
 * 	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *	@sa see also coupling::MultiMDMediator
 *  @author Niklas Wittmer
 */
template <class LinkedCell, unsigned int dim> class coupling::InstanceHandling {

public:
  /** Constructor:
   * 	@param mdConfig
   * 	@param mamicoConfig
   * 	@param multiMDService
   */
  InstanceHandling(simplemd::configurations::MolecularDynamicsConfiguration &mdConfig, coupling::configurations::MaMiCoConfiguration<dim> &mamicoConfig,
                   tarch::utils::MultiMDService<dim> &multiMDService)
      : _simpleMD(), _mdSolverInterface(), _mdConfig(mdConfig), _mamicoConfig(mamicoConfig), _multiMDService(multiMDService) {
    for (unsigned int i = 0; i < multiMDService.getLocalNumberOfMDSimulations(); i++) {
      _simpleMD.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(_mdConfig, _mamicoConfig
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                                                                            ,
                                                                                                            _multiMDService.getLocalCommunicator()
#endif
                                                                                                                ));

      if (_simpleMD[i] == nullptr) {
        std::cout << "ERROR InstanceHandling : _simpleMD [" << i << "] == NULL!" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      _simpleMD[i]->init(_multiMDService, _multiMDService.getGlobalNumberOfLocalMDSimulation(i));
    }
  }

  /** Destructor:
   */
  ~InstanceHandling() {
    for (unsigned int i = 0; i < _simpleMD.size(); ++i) {
      coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
      if (_simpleMD[i] != nullptr) {
        _simpleMD[i]->shutdown();
        delete _simpleMD[i];
        _simpleMD[i] = nullptr;
      }
      _mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().getMDSolverInterface();
    }
    _simpleMD.clear();
    for (auto &solverInterface : _mdSolverInterface) {
      if (solverInterface != nullptr) {
        delete solverInterface;
        solverInterface = nullptr;
      }
    }
    _mdSolverInterface.clear();
  }

  /** switches off the coupling between the new MD simulations and Macroscopic
   *solver and lets the MD simulations run t time steps starting from the time
   *step T to equilibrate. It should be called before switchOnCoupling()
   * 	@param T
   * 	@param t
   */
  void equilibrate(const unsigned int &t, const unsigned int &T) {
    for (auto &md : _simpleMD) {
      md->switchOffCoupling();
      md->simulateTimesteps(t, T);
    }
  }

  /** returns the vector of MD simulations
   *  @return  _simpleMD
   */
  auto &getSimpleMD() const { return _simpleMD; }

  /** Allocates Coupling interfaces
   * 	This method has to be called after switchOnCoupling()
   */
  void setMDSolverInterface() {

    for (unsigned int i = 0; i < _simpleMD.size(); ++i) {
      _mdSolverInterface.push_back(
          coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(_mdConfig, _mamicoConfig, _simpleMD[i]));
      if (_mdSolverInterface[i] == NULL) {
        std::cout << "ERROR InstanceHandling: mdSolverInterface[" << i << "] == NULL!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  /** Returns the MD solver onterface
   * 	@return _mdSolverInterface
   */
  auto &getMDSolverInterface() const { return _mdSolverInterface; }

  /** switches on the coupling between ALL new MD simulations and Macroscopic
   *solver after the new MD instances are equilibrated. It should be called
   *right after equilibrate() and before setMDSolverInterface()
   */
  void switchOnCoupling() {
    for (auto &simpleMD : _simpleMD) {
      simpleMD->switchOnCoupling();
    }
  }

  /** switches on the coupling between a specific new MD simulation with the
   *index number i and Macroscopic solver. It should be called right after
   *equilibrate() and before setMDSolverInterface()
   * 	@param i
   */
  void switchOnCoupling(const unsigned int &i) { _simpleMD[i]->switchOnCoupling(); }

  /** switches off the coupling between All new MD simulations and Macroscopic
   * solver.
   */
  void switchOffCoupling() {
    for (auto &simpleMD : _simpleMD) {
      simpleMD->switchOffCoupling();
    }
  }

  /** switches off the coupling between a new MD simulation with the index i and
   * Macroscopic solver.
   * 	@param i index number of the MD simulation
   */
  void switchOffCoupling(const unsigned int &i) { _simpleMD[i]->switchOffCoupling(); }

  /** Simulates t timesteps starting from current timestep T on all instances.
   * 	@param T
   * 	@param t
   */
  void simulateTimesteps(const unsigned int &t, unsigned int &T) {
    for (auto &simpleMD : _simpleMD) {
      simpleMD->simulateTimesteps(t, T);
    }
  }

  /** Saves the simulation result of the first MD instance as check point in the
   * file filestem
   * 	@param filestem
   * 	@param T
   */
  void writeCheckpoint(const std::string &filestem, const unsigned int &T) const {
    if (_simpleMD.size() > 0 && _simpleMD[0] != nullptr) {
      _simpleMD[0]->writeCheckpoint(filestem, T);
    }
  }

  /** Simulates t timesteps starting from current timestep T on all instances
   * 	and additionally uses MamicoInterfaceProvider for interfacing MD to FD.
   * 	@param t
   * 	@param T
   * 	@param multiMDCellService
   */
  void simulateTimesteps(const unsigned int &t, unsigned int &T, coupling::services::MultiMDCellService<LinkedCell, dim> &multiMDCellService) {
    for (unsigned int i = 0; i < _simpleMD.size(); ++i) {
      coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().setMacroscopicCellService(&multiMDCellService.getMacroscopicCellService(i));
      coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);

      if (_simpleMD[i] != nullptr) {
        _simpleMD[i]->simulateTimesteps(t, T);
      }
    }
  }

  /** Simulates t timesteps starting from current timestep T on all instances
   * but only performs simulation on one particular instance
   * 	@param t
   * 	@param T
   * 	@param i
   */
  void simulateTimesteps(const unsigned int &t, unsigned int &T, const unsigned int &i) { _simpleMD[i]->simulateTimesteps(t, T); }

  /** add a nullptr to the MD simulation vector and the vector of the MD solver
   * interface.
   */
  void addSimulationBlock() {
    _simpleMD.push_back(nullptr);
    _mdSolverInterface.push_back(nullptr);
  }

  /** rempve the last element of the MD simulation vector and the last element
   * of the vector of the MD solver interface.
   */
  void rmSimulationBlock() {
    _simpleMD.pop_back();
    _mdSolverInterface.pop_back();
  }

  /** adds one MS instance with the identifier localIndex to the slot "slot" and
   * return the corresponding MD solver interface. It initializes first an
   * instance of the slot "slot". Then set the _simpleMD[localIndex] to this new
   * instance.
   * 	@param slot
   * 	@param localIndex
   * 	@return _mdSolverInterface[localIndex]
   */
  coupling::interface::MDSolverInterface<LinkedCell, dim> *addMDSimulation(unsigned int slot, unsigned int localIndex) {
    auto *mdSim = coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(_mdConfig, _mamicoConfig
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                                                                    ,
                                                                                                    _multiMDService.getLocalCommunicator()
#endif
    );
    if (mdSim == NULL) {
      std::cout << "ERROR! coupling::InstanceHandling::addMDSimulation(): "
                   "mdSim == NULL!"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    mdSim->init(_multiMDService, slot);

    _simpleMD[localIndex] = mdSim;

    _mdSolverInterface[localIndex] =
        coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSolverInterface(_mdConfig, _mamicoConfig, _simpleMD[localIndex]);

    return _mdSolverInterface[localIndex];
  }

  /** removes one MS instance with the identifier index and delete the
   * corresponding MD solver interface.
   * 	@param index
   */
  void
  rmMDSimulation(const unsigned int &index) {
    if (_simpleMD[index] != nullptr) {
      _simpleMD[index]->shutdown();
      delete _simpleMD[index];
      _simpleMD[index] = nullptr;
    } else {
      std::cout << "WARNING coupling::InstanceHandling::rmMDSimulation() : "
                   "_simpleMD at index "
                << index << " == null!" << std::endl;
    }
    //_simpleMD.erase(_simpleMD.begin()+index);

    if (_mdSolverInterface[index] != nullptr) {
      delete _mdSolverInterface[index];
      _mdSolverInterface[index] = nullptr;
    } else {
      std::cout << "WARNING coupling::InstanceHandling::rmMDSimulation() : "
                   "_mdSolverInterface at index "
                << index << " == null!" << std::endl;
    }
    //_mdSolverInterface.erase(_mdSolverInterface.begin()+index);
  }

  /** sets single cell services in each MD simulation after initialising
   * macroscopic cell service for multi-MD case
   * 	@param multiMDCellService
   */
  void setMacroscopicCellServices(coupling::services::MultiMDCellService<LinkedCell, dim> &multiMDCellService) {
    for (unsigned int i = 0; i < _simpleMD.size(); ++i) {
      _simpleMD[i]->setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
    }
  }

private:
  std::vector<coupling::interface::MDSimulation *> _simpleMD;
  std::vector<coupling::interface::MDSolverInterface<LinkedCell, dim> *> _mdSolverInterface;

  simplemd::configurations::MolecularDynamicsConfiguration &_mdConfig;
  coupling::configurations::MaMiCoConfiguration<dim> &_mamicoConfig;

  const tarch::utils::MultiMDService<dim> &_multiMDService;
};

#endif //_INSTANCE_HANDLING_H_