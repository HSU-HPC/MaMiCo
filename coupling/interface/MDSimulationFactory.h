// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MDSIMULATIONFACTORY_H_
#define _MDSIMULATIONFACTORY_H_

#include <ctime>
#include <math.h>
#include <string>

// hacked: "currently, only 3D is supported!" (Defined before including headers in ./MDSimulationFactory/ using it)
#define MDSIMULATIONFACTORY_DIMENSION 3

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/configurations/MaMiCoConfiguration.h"
#include "coupling/interface/MDSimulation.h"
#include "coupling/interface/MamicoInterfaceProvider.h"
#include "simplemd/configurations/MolecularDynamicsConfiguration.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/utils/MultiMDService.h"

#if defined(SIMPLE_MD)
#include "coupling/interface/impl/SimpleMD/SimpleMDSolverInterface.h"
#include "simplemd/LinkedCell.h"
#define MY_LINKEDCELL simplemd::LinkedCell
#include "coupling/interface/impl/SimpleMD/SimpleMDSimulation.h"
#elif defined(LAMMPS_MD) || defined(LAMMPS_DPD)
#if defined(LAMMPS_MD)
#include "coupling/interface/impl/LAMMPS/LammpsMDSimulation.h"
#else // LAMMPS_DPD
#include "coupling/interface/impl/LAMMPS/LammpsDPDSimulation.h"
#endif
#define MY_LINKEDCELL LAMMPS_NS::MamicoCell
#elif defined(LS1_MARDYN)
#include "coupling/interface/impl/ls1/LS1MDSimulation.h"
#include "coupling/interface/impl/ls1/LS1MDSolverInterface.h"
#include "coupling/interface/impl/ls1/LS1RegionWrapper.h"
#include "coupling/interface/impl/ls1/LS1StaticCommData.h"
#define MY_LINKEDCELL ls1::LS1RegionWrapper
#endif

/** interface for different MD solvers.
 *  @author Philipp Neumann
 */
namespace coupling {
namespace interface {

/** factory to produced md simulation, md solver interface (for mamico) and the
 *coupling cell service using singleton pattern.
 *	@brief factory to produced md simulation, md solver interface (for
 *mamico) and the coupling cell service
 *  @author Philipp Neumann
 */
class SimulationAndInterfaceFactory {
public:
  /** @returns the SimulationAndInterfaceFactory object
   *	@note singleton pattern
   */
  static SimulationAndInterfaceFactory& getInstance() {
    static SimulationAndInterfaceFactory singleton;
    return singleton;
  }

  /** returns a pointer to the md simulation.
   *  @param configuration
   *  @param mamicoConfiguration
   *  @param localComm
   *	@remark This will create a new simulation which needs to be deleted at
   *the end
   */
  coupling::interface::MDSimulation* getMDSimulation(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                                                     const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                     ,
                                                     MPI_Comm localComm
#endif
  ) {
#if defined(SIMPLE_MD)
    return new coupling::interface::SimpleMDSimulation(configuration);
#elif defined(LAMMPS_MD)
    return new coupling::interface::LammpsMDSimulation(configuration, mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                       ,
                                                       localComm
#endif
    );
#elif defined(LAMMPS_DPD)
    return new coupling::interface::LammpsDPDSimulation(configuration, mamicoConfiguration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                        ,
                                                        localComm
#endif
    );
#elif defined(LS1_MARDYN)
    return new coupling::interface::LS1MDSimulation(configuration
#if (COUPLING_MD_PARALLEL == COUPLING_MD_YES)
                                                    ,
                                                    localComm
#endif
    );
#else
    std::cout << "ERROR MDSimulationFactory::getMDSimulation(): Unknown MD "
                 "simulation!"
              << std::endl;
    exit(EXIT_FAILURE);
    return NULL;
#endif
  }

  /** returns the MD solver interface. This method should be called AFTER
   * initialising the MD simulation AND AFTER the equilibration of MD. We thus
   * expect that getMDSimulationInterface() and switchOnCoupling() of the
   * MDSimulationInterface() have been called before.
   *  @param configuration
   *  @param mamicoConfiguration
   *  @param mdSimulation
   */
  coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>*
  getMDSolverInterface(const simplemd::configurations::MolecularDynamicsConfiguration& configuration,
                       const coupling::configurations::MaMiCoConfiguration<MDSIMULATIONFACTORY_DIMENSION>& mamicoConfiguration,
                       coupling::interface::MDSimulation* mdSimulation) {
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* mdSolverInterface = NULL;
#if defined(SIMPLE_MD)
    // for the simple MD code, we create a new md solver interface and also add
    // it to the mamico interface provider (the latter is not really required,
    // but makes the simulation state more consistent in the overall simulation)
    coupling::interface::SimpleMDSimulation* simpleMDSimulation = (coupling::interface::SimpleMDSimulation*)mdSimulation;
    if (simpleMDSimulation == NULL) {
      std::cout << "ERROR MDSimulationFactory::getMDSolverInterface(): Could "
                   "not cast to SimpleMDSimulation!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    mdSolverInterface = new coupling::interface::SimpleMDSolverInterface(
        simpleMDSimulation->getBoundaryTreatment(), simpleMDSimulation->getParallelTopologyService(), simpleMDSimulation->getMoleculeService(),
        simpleMDSimulation->getLinkedCellService(), simpleMDSimulation->getMolecularPropertiesService(),
        (simpleMDSimulation->getParallelTopologyService()).getLocalBoundaryInformation(), configuration.getSimulationConfiguration().getDt());
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(mdSolverInterface);
#elif defined(LAMMPS_MD)
    // as switchOnCoupling() should have been called before this method, the
    // fix-mamico should be already initialised. hence the MDSolverInterface of
    // the mamico interface provider should be initialised at this stage
    mdSolverInterface = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
#elif defined(LAMMPS_DPD)
    mdSolverInterface = coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
#elif defined(LS1_MARDYN)
    mdSolverInterface = new coupling::interface::LS1MDSolverInterface(mamicoConfiguration.getCouplingCellConfiguration().getCouplingCellSize(),
                                                                      mamicoConfiguration.getCouplingCellConfiguration().getNumberLinkedCellsPerCouplingCell());
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(mdSolverInterface);
#endif

    if (mdSolverInterface == NULL) {
      std::cout << "ERROR MDSimulationFactory::getMDSolverInterface(): "
                   "mdSolverInterface==NULL!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return mdSolverInterface;
  }

  /** shuts down the MD solver interfaces and deletes the interface if required
     (depending on the respective MD simulation) */
  void shutdownMDSolverInterface() {
#if defined(SIMPLE_MD)
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* interface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
    if (interface != NULL) {
      delete interface;
    }
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(NULL);
#elif defined(LAMMPS_MD)
// nop, since the fix-mamico takes care of deleting the MD solver interface
#elif defined(LAMMPS_DPD)
// nop
#elif defined(LS1_MARDYN)
    coupling::interface::MDSolverInterface<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>* interface =
        coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().getMDSolverInterface();
    if (interface != NULL) {
      delete interface;
    }
    coupling::interface::MamicoInterfaceProvider<MY_LINKEDCELL, MDSIMULATIONFACTORY_DIMENSION>::getInstance().setMDSolverInterface(NULL);
#endif
  }

private:
  /** Private constructor
   *  @note singelton pattern
   */
  SimulationAndInterfaceFactory() {}
  /** Private destructor
   *  @note singelton pattern
   */
  ~SimulationAndInterfaceFactory() {}
};

} // namespace interface
} // namespace coupling
#endif // _MDSIMULATIONFACTORY_H_
