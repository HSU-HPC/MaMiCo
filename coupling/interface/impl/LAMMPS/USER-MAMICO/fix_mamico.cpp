// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "fix_mamico.h"

LAMMPS_NS::FixMamico::FixMamico(LAMMPS* lmp, int argc, char** argv)
    : Fix(lmp, argc, argv), _lmp(lmp), _use2D(lmp->domain->dimension == 2), _timestepCounter(0) {
  if (argc != 6) {
    error->all(FLERR, "Specify the mamico fix as 'fix id group-id mamico maxCells seed "
                      "cutoffRadius' where maxCells is the maximum process-local number of "
                      "coupling cells (incl. ghost layer), seed is a random seed and "
                      "cutoffRadius is the cut-off radius in the respective LJ simulation");
  }

  int maxCells = atoi(argv[3]);
  int seed = atoi(argv[4]);
  double cutoffRadius = atof(argv[5]);
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "Max. number of local coupling cells to be allocated: " << atoi(argv[3]) << ", random seed: " << seed << std::endl;
#endif

  // access MaMiCo via InterfaceProvider
  if (_use2D) {
    MamicoLammpsMDSolverInterface<2>* mdSolverInterface2D = new MamicoLammpsMDSolverInterface<2>(_lmp, seed, maxCells, cutoffRadius);
    if (mdSolverInterface2D == NULL)
      error->all(FLERR, "FixMamico: mdSolverInterface2D==NULL!");
    coupling::interface::MamicoInterfaceProvider<MamicoCell, 2>::getInstance().setMDSolverInterface(mdSolverInterface2D);
  } else {
    MamicoLammpsMDSolverInterface<3>* mdSolverInterface3D = new MamicoLammpsMDSolverInterface<3>(_lmp, seed, maxCells, cutoffRadius);
    if (mdSolverInterface3D == NULL)
      error->all(FLERR, "FixMamico: mdSolverInterface3D==NULL!");
    coupling::interface::MamicoInterfaceProvider<MamicoCell, 3>::getInstance().setMDSolverInterface(mdSolverInterface3D);
  }
}

LAMMPS_NS::FixMamico::~FixMamico() {
  // destroy interface
  if (_use2D) {
    coupling::interface::MDSolverInterface<MamicoCell, 2>* interface =
        coupling::interface::MamicoInterfaceProvider<MamicoCell, 2>::getInstance().getMDSolverInterface();
    if (interface != NULL)
      delete interface;
    coupling::interface::MamicoInterfaceProvider<MamicoCell, 2>::getInstance().setMDSolverInterface(NULL);
  } else {
    coupling::interface::MDSolverInterface<MamicoCell, 3>* interface =
        coupling::interface::MamicoInterfaceProvider<MamicoCell, 3>::getInstance().getMDSolverInterface();
    if (interface != NULL)
      delete interface;
    coupling::interface::MamicoInterfaceProvider<MamicoCell, 3>::getInstance().setMDSolverInterface(NULL);
  }
}

void LAMMPS_NS::FixMamico::pre_force(int vflag) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "Enter FixMamico::pre_force()" << std::endl;
#endif
  // branching for 2D/3D
  if (_use2D) {
    modifyMDSystem<2>();
  } else {
    modifyMDSystem<3>();
  }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "Leave FixMamico::pre_force()" << std::endl;
#endif
}

void LAMMPS_NS::FixMamico::post_force(int vflag) {
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "Enter FixMamico::post_force()" << std::endl;
#endif
  // branching 2D/3D
  if (_use2D) {
    modifyMomentumAndTemperature<2>();
  } else {
    modifyMomentumAndTemperature<3>();
  }
#if (COUPLING_MD_DEBUG == COUPLING_MD_YES)
  std::cout << "Leave FixMamico::post_force()" << std::endl;
#endif
}

int LAMMPS_NS::FixMamico::setmask() {
  int mask = 0;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::POST_FORCE;
  return mask;
}
