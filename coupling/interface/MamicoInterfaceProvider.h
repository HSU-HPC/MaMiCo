// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_

#include "coupling/interface/MDSolverInterface.h"
#include "coupling/interface/MacroscopicSolverInterface.h"
#include "coupling/services/CouplingCellService.h"

namespace coupling {
namespace interface {

/** This is a singleton which returns and stores the interface implementations.
 *  It can be used if there is no other way to share access to a particular
 *interface object between MaMiCo and the MD simulation
 *	@brief a singleton which returns and stores the interface
 *implementations.
 *	@tparam LinkedCell cell type
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <class LinkedCell, int dim> class MamicoInterfaceProvider {
public:
  /** returns the MamicoInterfaceProvider object
   */
  static MamicoInterfaceProvider& getInstance() {
    static MamicoInterfaceProvider singleton;
    return singleton;
  }

  /** sets macroscopic solver interface
   *  @param macroscopicSolverInterface
   */
  void setMacroscopicSolverInterface(coupling::interface::MacroscopicSolverInterface<dim>* macroscopicSolverInterface) {
    _macroscopicSolverInterface = macroscopicSolverInterface;
  }

  /** returns acroscopic solver interface
   *  @return _macroscopicSolverInterface
   */
  coupling::interface::MacroscopicSolverInterface<dim>* getMacroscopicSolverInterface() { return _macroscopicSolverInterface; }

  /** sets MD solver interface
   *  @param mdSolverInterface
   */
  void setMDSolverInterface(coupling::interface::MDSolverInterface<LinkedCell, dim>* mdSolverInterface) { _mdSolverInterface = mdSolverInterface; }

  /** returns MD solver interface
   *  @return _mdSolverInterface
   */
  coupling::interface::MDSolverInterface<LinkedCell, dim>* getMDSolverInterface() { return _mdSolverInterface; }

  /** sets coupling cell service
   *  @return couplingCellService
   */
  void setCouplingCellService(coupling::services::CouplingCellService<dim>* couplingCellService) { _couplingCellService = couplingCellService; }

  /** returns coupling cell service
   *  @return _couplingCellService
   */
  coupling::services::CouplingCellService<dim>* getCouplingCellService() { return _couplingCellService; }

private:
  /** Private constructor, creation throgh a pointer and set functions
   *  @note singelton pattern
   */
  MamicoInterfaceProvider() : _macroscopicSolverInterface(NULL), _mdSolverInterface(NULL), _couplingCellService(NULL) {}
  /** Private destructor
   *  @note singelton pattern
   */
  ~MamicoInterfaceProvider() {
    _macroscopicSolverInterface = NULL;
    _mdSolverInterface = NULL;
    _couplingCellService = NULL;
  }

  coupling::interface::MacroscopicSolverInterface<dim>* _macroscopicSolverInterface;
  coupling::interface::MDSolverInterface<LinkedCell, dim>* _mdSolverInterface;
  coupling::services::CouplingCellService<dim>* _couplingCellService;
};

} // namespace interface
} // namespace coupling
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MAMICOINTERFACEPROVIDER_H_
