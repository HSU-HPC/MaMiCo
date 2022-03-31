// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVER_H_

#include "coupling/interface/impl/macroscopictestsolvers/TestMacroscopicSolver.h"

namespace coupling {
namespace interface { class VoidMacroscopicSolver; }
}

/** empty solver, does not do anything.
 * @author Philipp Neumann
 */
class coupling::interface::VoidMacroscopicSolver
    : public coupling::interface::TestMacroscopicSolver {
public:
  VoidMacroscopicSolver() : TestMacroscopicSolver() {}
  virtual ~VoidMacroscopicSolver() {}
  virtual void solve() {}
};

#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_VOIDMACROSCOPICSOLVER_H_
