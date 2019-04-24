#ifndef COMPRESSIBLENAVIERSTOKES_DENSITYCURRENT_H
#define COMPRESSIBLENAVIERSTOKES_DENSITYCURRENT_H

#include "Scenario.h"

namespace NavierStokes {

// TODO(Lukas) Refactor
class DensityCurrent : public Scenario {
 public:
  void initialValues(const double *const x, const PDE &ns,
                     Variables &vars) override;

  void source(const tarch::la::Vector<DIMENSIONS, double> &x, double t,
              const PDE &ns, const double *const Q, double *S) override;

  // Constants for dry air.
  // TODO(Lukas) Refactor these constants!
  // TODO(Lukas) Make sure that constants agree with each other.
  const double gamma = 1.4;
  const double Pr = 0.71;
  const double gasConstant = 287.058;
  const double c_p = 1.005 * 1000;
  const double c_v = 1 / (gamma - 1) * gasConstant;
  const double referencePressure = 10000;

  double getGamma() const override;

  double getPr() const override;

  double getC_v() const override;

  double getC_p() const override;

  double getGasConstant() const override;

  double getReferencePressure() const override;

  double getGravity() const override;

  double getBackgroundPotentialTemperature() const override;

  BoundaryType getBoundaryType(int faceId) override;

};

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_DENSITYCURRENT_H
