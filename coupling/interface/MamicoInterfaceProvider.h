// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_

#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/services/MacroscopicCellService.h"

namespace coupling {
namespace interface {

/** This is a singleton which returns and stores the interface implementations.
 *  It can be used if there is no other way to share access to a particular interface object between MaMiCo and the MD simulation.
 *  @author Philipp Neumann
 */
template<class LinkedCell,int dim>
class MamicoInterfaceProvider {
  public:
    static MamicoInterfaceProvider& getInstance(){
      static MamicoInterfaceProvider singleton; return singleton;
    }

    void setMacroscopicSolverInterface(coupling::interface::MacroscopicSolverInterface<dim> *macroscopicSolverInterface){
      _macroscopicSolverInterface = macroscopicSolverInterface;
    }
    coupling::interface::MacroscopicSolverInterface<dim>* getMacroscopicSolverInterface(){
      return _macroscopicSolverInterface;
    }
    void setMDSolverInterface(coupling::interface::MDSolverInterface<LinkedCell,dim>* mdSolverInterface){
      _mdSolverInterface = mdSolverInterface;
    }
    coupling::interface::MDSolverInterface<LinkedCell,dim>* getMDSolverInterface(){
      return _mdSolverInterface;
    }
    void setMacroscopicCellService(coupling::services::MacroscopicCellService<dim> *macroscopicCellService){
      _macroscopicCellService = macroscopicCellService;
    }
    coupling::services::MacroscopicCellService<dim>* getMacroscopicCellService(){
      return _macroscopicCellService;
    }

  private:
    MamicoInterfaceProvider(): _macroscopicSolverInterface(NULL),_mdSolverInterface(NULL),_macroscopicCellService(NULL){}
    ~MamicoInterfaceProvider(){ _macroscopicSolverInterface = NULL; _mdSolverInterface = NULL; _macroscopicCellService = NULL; }

    coupling::interface::MacroscopicSolverInterface<dim> *_macroscopicSolverInterface;
    coupling::interface::MDSolverInterface<LinkedCell,dim> *_mdSolverInterface;
    coupling::services::MacroscopicCellService<dim> *_macroscopicCellService;
};

}
}
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_

