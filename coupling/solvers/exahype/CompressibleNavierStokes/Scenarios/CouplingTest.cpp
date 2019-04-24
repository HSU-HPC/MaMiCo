#include "CouplingTest.h"
//#include "SmoothWave.h"

void NavierStokes::CouplingTest::initialValues(const double* const x,
                                             const PDE& ns, Variables& vars) {
  const auto initialZ = 1.0;
  TwoBubbles::initialValues(x, ns, vars, initialZ);
}
double NavierStokes::CouplingTest::getMolecularDiffusionCoeff() const {
  return 0.1; // smooth
}

double NavierStokes::CouplingTest::getQ0() const {
  return 0.0;
  //return 25.0;
}
void NavierStokes::CouplingTest::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  TwoBubbles::source(x, t, ns, Q, S);
  const auto k = 0.001; // decay rate
  ns.setZ(S, -k * ns.getZ(Q));
}

bool NavierStokes::CouplingTest::getUseAdvection() const {
  return true;
}