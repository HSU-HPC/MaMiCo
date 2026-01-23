// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

template <unsigned int dim> void LAMMPS_NS::FixMamico::modifyMDSystem() {
  coupling::services::CouplingCellService<dim>* couplingCellService =
      coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getCouplingCellService();
  const coupling::IndexConversion<dim>& indexConversion = couplingCellService->getIndexConversion();
  LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>* mdSolverInterface =
      (LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>*)coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMDSolverInterface();

// extract ghost atoms and update all mamico cell lists ------------------------
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  if (mdSolverInterface == NULL) {
    std::cout << "ERROR modifyMDSystem(): Could not cast MDSolverInterface!" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Sort particles into cells" << std::endl;
#endif
  mdSolverInterface->updateAllCells(indexConversion);

// call to coupling cell service functions ----------------------------------
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Process inner cells after MD timestep" << std::endl;
#endif
  couplingCellService->processInnerCouplingCellAfterMDTimestep();
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Distribute mass" << std::endl;
#endif
  couplingCellService->distributeMass(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico,Rank " << IDXS.getRank() << ": Plot every microscopic time step..." << std::endl;
#endif
  couplingCellService->plotEveryMicroscopicTimestep(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico,Rank " << IDXS.getRank() << ": Go to next time step..." << std::endl;
#endif
}

template <unsigned int dim> void LAMMPS_NS::FixMamico::modifyMomentumAndTemperature() {
  coupling::services::CouplingCellService<dim>* couplingCellService =
      coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getCouplingCellService();
  const coupling::IndexConversion<dim>& indexConversion = couplingCellService->getIndexConversion();
  LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>* mdSolverInterface =
      (LAMMPS_NS::MamicoLammpsMDSolverInterface<dim>*)coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMDSolverInterface();

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Apply thermostat" << std::endl;
#endif
  couplingCellService->applyTemperatureToMolecules(_timestepCounter);
  // update cells once again (sorting of molecules)
  mdSolverInterface->updateAllCells(indexConversion);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Modify momentum" << std::endl;
#endif
  couplingCellService->distributeMomentum(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << IDXS.getRank() << ": Apply boundary force" << std::endl;
#endif
  couplingCellService->applyBoundaryForce(_timestepCounter);

  // increment time step counter; since this method is called in post_force,
  // this is the right location for this purpose!
  _timestepCounter++;
}
