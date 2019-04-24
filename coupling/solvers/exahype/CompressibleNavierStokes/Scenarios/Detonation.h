#ifndef COMPRESSIBLENAVIERSTOKES_DETONATION_H
#define COMPRESSIBLENAVIERSTOKES_DETONATION_H

#include "Scenario.h"

namespace NavierStokes {

class Detonation : public Scenario {
public:
  void initialValues(const double *const x, const PDE &ns,
                     Variables &vars) final override;
  double getMolecularDiffusionCoeff() const override;
  double getQ0() const override;
  double getGamma() const override;
  double getPr() const override;
  double getC_v() const override;
  double getC_p() const override;
  double getGasConstant() const override;
  BoundaryType getBoundaryType(int faceId) override;
  bool getUseAdvection() const override;
  void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
              const PDE& ns, const double* const Q, double* S) override;
};
}
#endif //COMPRESSIBLENAVIERSTOKES_DETONATION_H
