#ifndef COMPRESSIBLENAVIERSTOKES_LIDDRIVENCAVITY_H
#define COMPRESSIBLENAVIERSTOKES_LIDDRIVENCAVITY_H

#include "Scenario.h"

namespace NavierStokes {
    class LidDrivenCavity : public Scenario {
        void initialValues(const double *const x, const PDE &ns,
                           Variables &vars) override;

        BoundaryType getBoundaryType(int faceId) override;
    };

}

#endif //COMPRESSIBLENAVIERSTOKES_LIDDRIVENCAVITY_H
