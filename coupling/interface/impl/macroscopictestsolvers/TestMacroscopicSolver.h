// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVER_H_

namespace coupling {
namespace interface { class TestMacroscopicSolver; }
}

/** interface for different test-macroscopic solvers.
 *  @author Philipp Neumann
 */
class coupling::interface::TestMacroscopicSolver {
public:
  TestMacroscopicSolver() {}
  virtual ~TestMacroscopicSolver() {}

  virtual void solve() = 0;
};
#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_TESTMACROSCOPICSOLVER_H_
