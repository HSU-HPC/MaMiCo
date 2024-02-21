#pragma once

#include "coupling/interface/MDSimulation.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#include "ls1/src/Domain.h"
#include "ls1/src/Simulation.h"
#include "ls1/src/plugins/MamicoCoupling.h"

namespace coupling {
namespace interface {
class LS1MDSimulation : public coupling::interface::MDSimulation {
private:
  const simplemd::configurations::MolecularDynamicsConfiguration& _configuration;
  Simulation* simulation;          // cannot name this _simulation, a global preprocessor marco with the name _simulation expands to *global_simulation
  MamicoCoupling* ls1MamicoPlugin; // the plugin is only initialized after the simulation object reads xml, so cannot use it before that point
  bool internalCouplingState;
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
  MPI_Comm comm;
#endif

public:
  LS1MDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                  ,
                  MPI_Comm localComm
#endif
                  )
      : coupling::interface::MDSimulation(), _configuration(configuration) {

#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
    comm = localComm;
    LS1StaticCommData::getInstance().setLocalCommunicator(comm); // needs to be done before creating the simulation object
    const tarch::la::Vector<3, unsigned int> numberProcs = _configuration.getMPIConfiguration().getNumberOfProcesses();
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(0, numberProcs[0]);
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(1, numberProcs[1]);
    LS1StaticCommData::getInstance().setDomainGridDecompAtDim(2, numberProcs[2]);
#endif

    simulation = new Simulation();
    global_simulation = simulation;
    simulation->disableFinalCheckpoint();
    internalCouplingState = false;
    ls1MamicoPlugin = nullptr;
  }
  virtual ~LS1MDSimulation() {
    if (simulation != nullptr) {
      delete simulation;
      simulation = nullptr;
    }
  }
  /** switches coupling on/off*/
  virtual void switchOffCoupling() override {
    // coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOff();
    internalCouplingState = false;
    simulation->getDomain()->thermostatOn();
    simulation->getDomain()->setExplosionHeuristics(true);
    if (ls1MamicoPlugin != nullptr)
      ls1MamicoPlugin->switchOffCoupling();
  }
  virtual void switchOnCoupling() override {
    // coupling::interface::LS1MamicoCouplingSwitch::getInstance().setCouplingStateOn();
    internalCouplingState = true;
    simulation->getDomain()->thermostatOff();
    simulation->getDomain()->setExplosionHeuristics(false);
    if (ls1MamicoPlugin != nullptr)
      ls1MamicoPlugin->switchOnCoupling();
  }

  /** simulates numberTimesteps time steps and starts at time step no.
   * firstTimestep*/
  virtual void simulateTimesteps(const unsigned int& numberTimesteps, const unsigned int& firstTimestep) override {
    global_simulation = simulation;
    for (unsigned int i = 0; i < numberTimesteps; i++) {
      simulation->simulateOneTimestep();
    }
  }
  /** simulates a single time step*/
  // virtual void simulateTimestep(const unsigned int &thisTimestep ){const
  // unsigned int steps=1; simulateTimesteps(thisTimestep,steps);} TODO BUG
  virtual void sortMoleculesIntoCells() override {}

  virtual void setCouplingCellService(coupling::services::MacroscopicCellService<MDSIMULATIONFACTORY_DIMENSION>* macroscopicCellService) override {
    // coupling::interface::MamicoInterfaceProvider<ls1::LS1RegionWrapper, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setCouplingCellService(
    //     macroscopicCellService);
    global_simulation = simulation;
    PluginBase* searchedPlugin = simulation->getPlugin("MamicoCoupling");
    if (searchedPlugin == nullptr) {
      std::cout << "ERROR: MaMiCo plugin not found!" << std::endl;
      exit(EXIT_FAILURE);
    }
    ls1MamicoPlugin = dynamic_cast<MamicoCoupling*>(searchedPlugin);
    if (ls1MamicoPlugin != nullptr) {
      ls1MamicoPlugin->setMamicoMacroscopicCellService(macroscopicCellService);
    } else {
      std::cout << "ERROR: Cast to Mamico plugin unsuccessful!" << std::endl;
      exit(EXIT_FAILURE);
    }
    // since this is the first time the plugin is accessed, set whatever preexisting coupling variable we had here for the first time
    if (internalCouplingState)
      ls1MamicoPlugin->switchOnCoupling();
    else
      ls1MamicoPlugin->switchOffCoupling();
  }
  virtual void init() override {
    global_simulation = simulation;
    // parse file
    const std::string filename = coupling::interface::LS1StaticCommData::getInstance().getConfigFilename();
    simulation->readConfigFile(filename);
    simulation->getDomain()->thermostatOff();
    simulation->getDomain()->setExplosionHeuristics(false);
    // after this point the mamico plugin exists and is accessible
    simulation->prepare_start();
    simulation->preSimLoopSteps();
  }
  virtual void init(const tarch::utils::MultiMDService<MDSIMULATIONFACTORY_DIMENSION>& multiMDService, unsigned int localMDSimulation) override { init(); }
  virtual void shutdown() override {
    global_simulation = simulation;
    simulation->markSimAsDone();
    simulation->postSimLoopSteps();
    simulation->finalize();
  }
  virtual void writeCheckpoint(const std::string& filestem, const unsigned int& t) override {
    // configure through ls1 config file, using plugins
  }
};
} // namespace interface
} // namespace coupling
