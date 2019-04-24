#ifndef COMPRESSIBLENAVIERSTOKES_ENTROPYWAVE_H
#define COMPRESSIBLENAVIERSTOKES_ENTROPYWAVE_H

#include "Scenario.h"

namespace NavierStokes {
class EntropyWave : public Scenario {
 public:
  void analyticalSolution(const double* const x, double t, const PDE& ns,
                          Variables& vars, double* gradState) override;
  BoundaryType getBoundaryType(int faceId) override;
};
}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_ENTROPYWAVE_H
