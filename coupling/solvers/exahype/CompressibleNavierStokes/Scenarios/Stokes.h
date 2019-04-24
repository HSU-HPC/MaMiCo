#ifndef NAVIERSTOKES_STOKES_H
#define NAVIERSTOKES_STOKES_H

#include "Scenario.h"

namespace NavierStokes {
class Stokes : public Scenario {
  void initialValues(const double* const x, const PDE& ns,
                     Variables& vars) final override;
  void analyticalSolution(const double* const x, const double t, const PDE& ns,
                          Variables& vars, double* gradState) final override;
  BoundaryType getBoundaryType(int faceId) final override;
};

}  // namespace NavierStokes

#endif  // NAVIERSTOKES_STOKES_H
