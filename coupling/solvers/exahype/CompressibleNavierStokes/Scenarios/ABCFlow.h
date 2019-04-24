#ifndef COMPRESSIBLENAVIERSTOKES_ABCFLOW_H
#define COMPRESSIBLENAVIERSTOKES_ABCFLOW_H

#include "Scenario.h"

namespace NavierStokes {

class ABCFlow : public Scenario {
  void analyticalSolution(const double* const x, const double t, const PDE& ns,
                          Variables& vars, double* gradState) final override;
  BoundaryType getBoundaryType(int faceId) final override;
};
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_ABCFLOW_H
