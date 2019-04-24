#ifndef NAVIERSTOKES_SODSHOCKTUBE_H
#define NAVIERSTOKES_SODSHOCKTUBE_H

#include "Scenario.h"

namespace NavierStokes {
class SodShockTube : public Scenario {
  void initialValues(const double* const x, const PDE& ns,
                     Variables& vars) final override;
};

}  // namespace NavierStokes

#endif  // NAVIERSTOKES_SODSHOCKTUBE_H
