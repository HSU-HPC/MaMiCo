// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

template <unsigned int dim> void LAMMPS_NS::FixMamico::modifyMDSystem() {
  coupling::services::MacroscopicCellService<dim> *macroscopicCellService =
      coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMacroscopicCellService();
  const coupling::IndexConversion<dim> &indexConversion = macroscopicCellService->getIndexConversion();
  LAMMPS_NS::MamicoLammpsMDSolverInterface<dim> *mdSolverInterface =
      (LAMMPS_NS::MamicoLammpsMDSolverInterface<dim> *)coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMDSolverInterface();

// extract ghost atoms and update all mamico cell lists ------------------------
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  if (mdSolverInterface == NULL) {
    std::cout << "ERROR modifyMDSystem(): Could not cast MDSolverInterface!" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Sort particles into cells" << std::endl;
#endif
  mdSolverInterface->updateAllCells(indexConversion);

// call to macroscopic cell service functions ----------------------------------
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Process inner cells after MD timestep" << std::endl;
#endif
  macroscopicCellService->processInnerCouplingCellAfterMDTimestep();
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Distribute mass" << std::endl;
#endif
  macroscopicCellService->distributeMass(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico,Rank " << indexConversion.getThisRank() << ": Plot every microscopic time step..." << std::endl;
#endif
  macroscopicCellService->plotEveryMicroscopicTimestep(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico,Rank " << indexConversion.getThisRank() << ": Go to next time step..." << std::endl;
#endif
}

template <unsigned int dim> void LAMMPS_NS::FixMamico::modifyMomentumAndTemperature() {
  coupling::services::MacroscopicCellService<dim> *macroscopicCellService =
      coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMacroscopicCellService();
  const coupling::IndexConversion<dim> &indexConversion = macroscopicCellService->getIndexConversion();
  LAMMPS_NS::MamicoLammpsMDSolverInterface<dim> *mdSolverInterface =
      (LAMMPS_NS::MamicoLammpsMDSolverInterface<dim> *)coupling::interface::MamicoInterfaceProvider<MamicoCell, dim>::getInstance().getMDSolverInterface();

#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Apply thermostat" << std::endl;
#endif
  macroscopicCellService->applyTemperatureToMolecules(_timestepCounter);
  // update cells once again (sorting of molecules)
  mdSolverInterface->updateAllCells(indexConversion);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Modify momentum" << std::endl;
#endif
  macroscopicCellService->distributeMomentum(_timestepCounter);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "FixMamico, Rank " << indexConversion.getThisRank() << ": Apply boundary force" << std::endl;
#endif
  macroscopicCellService->applyBoundaryForce(_timestepCounter);

  // increment time step counter; since this method is called in post_force,
  // this is the right location for this purpose!
  _timestepCounter++;
}
