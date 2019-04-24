#ifndef COMPRESSIBLENAVIERSTOKES_COUPLINGTEST_H
#define COMPRESSIBLENAVIERSTOKES_COUPLINGTEST_H

#include "Scenario.h"
#include "TwoBubbles.h"

namespace NavierStokes {
class CouplingTest : public TwoBubbles {
public:
  void initialValues(const double *const x, const PDE &ns,
                     Variables &vars) final override;
  double getMolecularDiffusionCoeff() const override;
  double getQ0() const override;
  bool getUseAdvection() const override;
  void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
              const PDE& ns, const double* const Q, double* S) override;
};
}  // namespace NavierStokes

#endif //COMPRESSIBLENAVIERSTOKES_COUPLINGTEST_H
