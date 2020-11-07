// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _INSTANCE_HANDLING_H_
#define _INSTANCE_HANDLING_H_

#include "tarch/utils/MultiMDService.h"
#include "coupling/services/MultiMDCellService.h"
#include "coupling/interface/MDSimulationFactory.h"

namespace coupling {
    class InstanceHandling;
}


class coupling::InstanceHandling {

public:
  InstanceHandling(
      simplemd::configurations::MolecularDynamicsConfiguration& mdConfig,
      coupling::configurations::MaMiCoConfiguration& mamicoConfig,
      const tarch::utils::MultiMDService& multiMDService
    ) : _simpleMD(), _mdSolverInterface(), 
        _mdConfig(mdConfig), _mamicoConfig(mamicoConfig), _multiMDService(multiMDService)
    {
      for(unsigned int i=0;i<multiMDService.getLocalNumberOfMDSimulations()) {
        _simpleMD.push_back(
            coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
              _mdConfig, _mamicoConfig
              #if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
              , _multiMDService.getLocalCommunicator()
              #endif
        ));

        if(_simpleMD[i]==nullptr) {
          std::cout << "ERROR InstanceHandling : _simpleMD [" << i << "] == NULL!" << std::endl;
          std::exit(EXIT_FAILURE);
        }

        _simpleMD[i]->init(_multiMDService,_multiMDService.getGlobalNumberOfLocalMDSimulation(i));

        // Allocate coupling interfaces
        _mdSolverInterface.push_back(
          coupling::interface::SimulationAndInterfaceFactory::getInstance()
            .getMDSolverInterface(_mdConfig, _mamicoConfig, _simpleMD[i])
        );
        if (_mdSolverInterface[i] == NULL){
          std::cout << "ERROR InstanceHandling: mdSolverInterface[" << i << "] == NULL!" << std::endl; 
          exit(EXIT_FAILURE);
        }
      } 
    }

    auto & getSimpleMD() const {
      return _simpleMD;
    }

    auto & getMDSolverInterface() const {
      return _mdSolverInterface;
    }

    void switchOnCoupling() {
      for(auto & simpleMD : _simpleMD) {
        simpleMD.switchOnCoupling();
      }
    }
    void switchOnCoupling(const unsigned int & i) {
      _simpleMD[i].switchOnCoupling();
    }

    void switchOffCoupling() {
        for(auto & simpleMD : _simpleMD) {
            simpleMD.switchOnCoupling();
        }
    }
    void switchOffCoupling(const unsigned int & i) {
        _simpleMD[i].switchOffCoupling();
    }

    // Simulates t timesteps starting frum current timestep T on all instances
    void simulateTimesteps(const unsigned int & t, unsigned int & T) {
        for(auto & simpleMD : _simpleMD) {
            simpleMD.simulateTimesteps(t, T);
        }
    }

    // Same as above
    // but additionally uses MamicoInterfaceProvider for interfacing MD to FD
    void simulateTimesteps(const unsigned int & t, unsigned int & T, 
                            const coupling::services::MultiMDCellService<MY_LINKEDCELL,3> & multiMDCellService) {
        for(unsigned int i =0;i<simplemd.size();++i) {
          coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>
            ::getInstance.setMacroscopicCellService(multiMDCellService.getMacroscopicCellService(i));
          coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,3>
            ::getInstance.setMDSolverInterface(_mdSolverInterface[i]);

          _simpleMD[i].simulateTimesteps(_mdConfig.getSimulationConfiguration().getNumberOfTimesteps());
        }
    }

    // Same as simulateTimesteps(t, T) but only performs simulation on one particular instance
    void simulateTimesteps(const unsigned int & t, unsigned int & T, const unsigned int & i) {
        _simpleMD[i].simulateTimesteps(t,T);
    }



private:
    std::vector<coupling::interface::MDSimulation&>& _simpleMD;
    std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, 3>& >& _mdSolverInterface;

    simplemd::configurations::MolecularDynamicsConfiguration& _mdConfig;
    coupling::configurations::MaMiCoConfiguration& _mamicoConfig;

    const tarch::utils::MultiMDService& _multiMDService;

};

#endif _INSTANCE_HANDLING_H_