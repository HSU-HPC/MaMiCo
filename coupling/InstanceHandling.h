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
    template<class LinkedCell, unsigned int dim> class InstanceHandling;
}

template<class LinkedCell, unsigned int dim>
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

        
      } 
    }

    ~InstanceHandling() {
      for(unsigned int i=0;i<_simpleMD.size();++i) {
        coupling::interface::MamicoInterfaceProvider<LinkedCell,dim>::getInstance().setMDSolverInterface(_mdSolverInterface[i]);
        if(_simpleMD[i] != nullptr) {
          _simpleMD[i]->shutdown();
          delete _simpleMD[i];
          _simpleMD[i] = nullptr;
        }
        _mdSolverInterface[i] = coupling::interface::MamicoInterfaceProvider<LinkedCell, dim>::getInstance().getMDSolverInterface();
      }
      _simpleMD.clear();
      for(auto & solverInterface : _mdSolverInterface) {
        if(solverInterface != nullptr) {
          delete solverInterface;
          solverInterface = nullptr;
        }
      }
      _mdSolverInterface.clear();

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

    /** Allocate Coupling interfaces
     * This method has to be called after switchOnCoupling()
     */
    void setMDSolverInterface() {
        
        for(unsigned int i=0;i<_simpleMD.size();++i) {
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

    void writeCheckpoint(const std::string & filestem, const unsigned int & T) const {
      if(_simpleMD.size() > 0 && _simpleMD[0] != nullptr) {
        _simpleMD[0]->writeCheckpoint(filestem, T);
      }
    }

    // Same as above
    // but additionally uses MamicoInterfaceProvider for interfacing MD to FD
    void simulateTimesteps(const unsigned int & t, unsigned int & T, 
                            coupling::services::MultiMDCellService<LinkedCell,dim> & multiMDCellService) {
        for(unsigned int i =0;i<_simpleMD.size();++i) {
          coupling::interface::MamicoInterfaceProvider<LinkedCell,dim>
            ::getInstance().setMacroscopicCellService(&multiMDCellService.getMacroscopicCellService(i));
          coupling::interface::MamicoInterfaceProvider<LinkedCell,dim>
            ::getInstance().setMDSolverInterface(_mdSolverInterface[i]);

          if(_simpleMD[i] != nullptr) {
            _simpleMD[i]->simulateTimesteps(t,T);
          }

        }
    }

    // Same as simulateTimesteps(t, T) but only performs simulation on one particular instance
    void simulateTimesteps(const unsigned int & t, unsigned int & T, const unsigned int & i) {
        _simpleMD[i]->simulateTimesteps(t,T);
    }

    void addSimulationBlock() {
      _simpleMD.push_back(nullptr);
      _mdSolverInterface.push_back(nullptr);
    }

    void rmSimulationBlock() {
      _simpleMD.pop_back();
      _mdSolverInterface.pop_back();
    }

    coupling::interface::MDSolverInterface<LinkedCell,dim>* addMDSimulation(unsigned int slot, unsigned int localIndex) {
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

      _simpleMD[localIndex] = mdSim;

      _mdSolverInterface[localIndex] = coupling::interface::SimulationAndInterfaceFactory::getInstance()
                                        .getMDSolverInterface(_mdConfig, _mamicoConfig, _simpleMD[localIndex]);

      return _mdSolverInterface[localIndex];
    }

    void rmMDSimulation(const unsigned int & index) {
      if(_simpleMD[index] != nullptr) {
        _simpleMD[index]->shutdown();
        delete _simpleMD[index];
        _simpleMD[index] = nullptr;
      } else {
        std::cout << "WARNING coupling::InstanceHandling::rmMDSimulation() : _simpleMD at index " << index << " == null!" << std::endl;
      }
      //_simpleMD.erase(_simpleMD.begin()+index);

      if(_mdSolverInterface[index] != nullptr) {
        delete _mdSolverInterface[index];
        _mdSolverInterface[index] = nullptr;
      } else {
        std::cout << "WARNING coupling::InstanceHandling::rmMDSimulation() : _mdSolverInterface at index " << index << " == null!" << std::endl;
      }
      //_mdSolverInterface.erase(_mdSolverInterface.begin()+index);
    }

    void setMacroscopicCellServices(coupling::services::MultiMDCellService<LinkedCell, dim> & multiMDCellService) {
      for(unsigned int i=0;i<_simpleMD.size();++i) {
        _simpleMD[i]->setMacroscopicCellService(&(multiMDCellService.getMacroscopicCellService(i)));
      }
    }

private:
    std::vector<coupling::interface::MDSimulation*> _simpleMD;
    std::vector<coupling::interface::MDSolverInterface<LinkedCell, dim>* > _mdSolverInterface;

    simplemd::configurations::MolecularDynamicsConfiguration& _mdConfig;
    coupling::configurations::MaMiCoConfiguration<dim>& _mamicoConfig;

    const tarch::utils::MultiMDService<dim>& _multiMDService;
    

};

#endif //_INSTANCE_HANDLING_H_