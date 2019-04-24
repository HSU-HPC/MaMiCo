#ifndef NAVIERSTOKES_SCENARIO_H
#define NAVIERSTOKES_SCENARIO_H

#include "../PDE.h"
#include "tarch/la/Vector.h"
#include <array>

namespace NavierStokes {
enum class BoundaryType {analytical, wall, freeSlipWall, hydrostaticWall, movingWall };

class Scenario {
 public:
  Scenario();
  // TODO(Lukas) Mark these all as const!
  virtual void initialValues(const double* const x, const PDE& ns,
                             Variables& vars);
  virtual void analyticalSolution(const double* const x, double t,
                                  const PDE& ns, Variables& vars,
                                  double* gradState);
  virtual void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
                      const PDE& ns, const double* const Q, double* S);
  virtual BoundaryType getBoundaryType(int faceId);

  virtual double getGamma() const;
  virtual double getPr() const;
  virtual double getC_v() const;
  virtual double getC_p() const;
  virtual double getGasConstant() const;
  virtual double getReferencePressure() const;
  virtual double getGravity() const;
  virtual double getBackgroundPotentialTemperature() const;
  virtual bool getUseAdvection() const;

  virtual double getMolecularDiffusionCoeff() const;
  virtual double getQ0() const;

  const double gamma = 1.4;
  const double Pr = 0.7;
  const double c_v = 1;
  const double c_p = c_v * gamma;
  const double gasConstant = c_p - c_v;
  const double referencePressure = 10000;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SCENARIO_H
