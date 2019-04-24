#ifndef NAVIERSTOKES_SMOOTHWAVE_H
#define NAVIERSTOKES_SMOOTHWAVE_H

#include "Scenario.h"

namespace NavierStokes {
class SmoothWave : public Scenario {
public:
  void initialValues(const double *const x, const PDE &ns,
                     Variables &vars) final override;
  BoundaryType getBoundaryType(int faceId) override;
};
}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SMOOTHWAVE_H
