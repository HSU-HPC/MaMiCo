#ifndef COMPRESSIBLENAVIERSTOKES_COUETTE_H
#define COMPRESSIBLENAVIERSTOKES_COUETTE_H

#include "Scenario.h"

namespace NavierStokes {
    class Couette : public Scenario {
        void initialValues(const double *const x, const PDE &ns,
                           Variables &vars) override;
        
        void analyticalSolution(const double* const x, double t, const PDE& ns,
                          Variables& vars, double* gradState) override;

        BoundaryType getBoundaryType(int faceId) override;
    private:
    	double getVelocity(double z);
    };

}

#endif //COMPRESSIBLENAVIERSTOKES_COUETTE_H
