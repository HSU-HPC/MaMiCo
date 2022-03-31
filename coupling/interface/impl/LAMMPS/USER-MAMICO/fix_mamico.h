// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifdef FIX_CLASS

FixStyle(mamico, FixMamico)

#else

#ifndef LMP_FIX_MAMICO_H
#define LMP_FIX_MAMICO_H

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "lammps.h"

#include <cstdlib>
#include <iostream>

#include "coupling/CouplingMDDefinitions.h"
#include "coupling/interface/MamicoInterfaceProvider.h"

#include "mamico_cell.h"
#include "mamico_lammps_md_solver_interface.h"
#include "sorting.h"

namespace LAMMPS_NS {
/** fix to hook in and add MaMiCo functionality.
 *  @author Philipp Neumann
 */
class FixMamico : public Fix {
public:
  FixMamico(LAMMPS *lmp, int argc, char **argv);
  ~FixMamico();

  void pre_force(int vflag);
  void post_force(int vflag);

  int setmask();

private:
  // carries out the coupling to MD system and inserts/deletes particles. Also
  // triggers plotting at MD scale.
  // We needed to templatize this part to reduce code duplication w.r.t. to
  // 2D/3D.
  // Called from pre_force() callback.
  template <unsigned int dim> void modifyMDSystem();

  // modifies momentum, applies thermostat and triggers incrementation of time
  // step counter.
  // Called from post_force() callback.
  template <unsigned int dim> void modifyMomentumAndTemperature();

  template <unsigned int dim> void sortAtomPositionsIntoCells(const coupling::IndexConversion<dim> &indexConversion);

  LAMMPS *_lmp;                  // ptr to lammps
  bool _use2D;                   // true, if this is a 2D simulation
  unsigned int _timestepCounter; // counts the total number of time steps where
                                 // this fix is applied
};
} // namespace LAMMPS_NS

#include "fix_mamico_template_functions.h"

#endif // LMP_FIX_MAMICO_H
#endif
