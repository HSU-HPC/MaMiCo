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
    template<unsigned int dim> class InstanceHandling;
}

template<unsigned int dim>
class coupling::InstanceHandling {

public:
  InstanceHandling(
      simplemd::configurations::MolecularDynamicsConfiguration& mdConfig,
      coupling::configurations::MaMiCoConfiguration<dim>& mamicoConfig,
      tarch::utils::MultiMDService<dim>& multiMDService
    ) : _simpleMD(), _mdSolverInterface(), 
        _mdConfig(mdConfig), _mamicoConfig(mamicoConfig), _multiMDService(multiMDService)
    {
      for(unsigned int i=0;i<multiMDService.getLocalNumberOfMDSimulations();i++) {
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

    void equilibrate(const unsigned int & t, const unsigned int & T) {
      for(auto & md : _simpleMD) {
        md->switchOffCoupling();
        md->simulateTimesteps(t, T);  
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
        simpleMD->switchOnCoupling();
      }
    }
    void switchOnCoupling(const unsigned int & i) {
      _simpleMD[i]->switchOnCoupling();
    }

    void switchOffCoupling() {
        for(auto & simpleMD : _simpleMD) {
            simpleMD->switchOnCoupling();
        }
    }
    void switchOffCoupling(const unsigned int & i) {
        _simpleMD[i]->switchOffCoupling();
    }

    // Simulates t timesteps starting frum current timestep T on all instances
    void simulateTimesteps(const unsigned int & t, unsigned int & T) {
        for(auto & simpleMD : _simpleMD) {
            simpleMD->simulateTimesteps(t, T);
        }
    }

    // Same as above
    // but additionally uses MamicoInterfaceProvider for interfacing MD to FD
    void simulateTimesteps(const unsigned int & t, unsigned int & T, 
                            coupling::services::MultiMDCellService<MY_LINKEDCELL,dim> & multiMDCellService) {
        for(unsigned int i =0;i<_simpleMD.size();++i) {
          coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,dim>
            ::getInstance().setMacroscopicCellService(&multiMDCellService.getMacroscopicCellService(i));
          coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL,dim>
            ::getInstance().setMDSolverInterface(_mdSolverInterface[i]);

          _simpleMD[i]->simulateTimesteps(t,T);
        }
    }

    // Same as simulateTimesteps(t, T) but only performs simulation on one particular instance
    void simulateTimesteps(const unsigned int & t, unsigned int & T, const unsigned int & i) {
        _simpleMD[i]->simulateTimesteps(t,T);
    }

    coupling::interface::MDSolverInterface<MY_LINKEDCELL,dim>* addMDSimulation(unsigned int slot) {
      auto * mdSim = coupling::interface::SimulationAndInterfaceFactory::getInstance().getMDSimulation(
        _mdConfig, _mamicoConfig
        #if (COUPLING_MD_PARALLEL==COUPLING_MD_YES)
        , _multiMDService.getLocalCommunicator()
        #endif
      );
      if(mdSim == NULL) {
          std::cout << "ERROR! coupling::InstanceHandling::addMDSimulation(): mdSim == NULL!" << std::endl;
          std::exit(EXIT_FAILURE);
      }

      mdSim->init(_multiMDService, slot);

      _simpleMD.push_back(mdSim);

      _mdSolverInterface.push_back(coupling::interface::SimulationAndInterfaceFactory::getInstance()
                                    .getMDSolverInterface(_mdConfig, _mamicoConfig, _simpleMD[_simpleMD.size()-1]));

      return _mdSolverInterface[_mdSolverInterface.size()-1];
    }

    void rmMDSimulation(const unsigned int & index) {
      _simpleMD[index]->shutdown();
      delete _simpleMD[index];
      _simpleMD[index] = nullptr;
      _simpleMD.erase(_simpleMD.begin()+index);

      delete _mdSolverInterface[index];
      _mdSolverInterface[index] = nullptr;
      _mdSolverInterface.erase(_mdSolverInterface.begin()+index);
    }

private:
    std::vector<coupling::interface::MDSimulation*> _simpleMD;
    std::vector<coupling::interface::MDSolverInterface<MY_LINKEDCELL, dim>* > _mdSolverInterface;

    simplemd::configurations::MolecularDynamicsConfiguration& _mdConfig;
    coupling::configurations::MaMiCoConfiguration<dim>& _mamicoConfig;

    const tarch::utils::MultiMDService<dim>& _multiMDService;

};

#endif //_INSTANCE_HANDLING_H_