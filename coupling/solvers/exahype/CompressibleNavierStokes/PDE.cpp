#include "PDE.h"
#include "Scenarios/Scenario.h"

NavierStokes::PDE::PDE() :
  PDE(0.1, 10000, 1.4, 0.7, 1, 1.4, 0.4) {
}

NavierStokes::PDE::PDE(double referenceViscosity, NavierStokes::Scenario &scenario, bool useGravity,
        bool useBackgroundState) :
  referenceViscosity(referenceViscosity),
  referencePressure(scenario.getReferencePressure()),
  gamma(scenario.getGamma()),
  Pr(scenario.getPr()),
  c_v(scenario.getC_v()),
  c_p(scenario.getC_p()),
  gasConstant(scenario.getGasConstant()),
  q0(scenario.getQ0()),
  molecularDiffusionCoeff(scenario.getMolecularDiffusionCoeff()),
  useGravity(useGravity),
  useBackgroundState(useBackgroundState),
  gravitation(scenario.getGravity()),
  useAdvection(scenario.getUseAdvection()) {
    assertion3(q0 == 0.0 || useAdvection, q0, molecularDiffusionCoeff, useAdvection);
}

NavierStokes::PDE::PDE(double referenceViscosity, double referencePressure, double gamma, double Pr,
        double c_v, double c_p, double gasConstant) :
  referenceViscosity(referenceViscosity),
  referencePressure(referencePressure),
  gamma(gamma),
  Pr(Pr),
  c_v(c_v),
  c_p(c_p),
  gasConstant(gasConstant) {
    q0 = 0.0;
    molecularDiffusionCoeff = 0.0;
    useAdvection = false;
    useGravity = false;
    useBackgroundState = false;
    gravitation = 0.0;
}

double NavierStokes::PDE::getZ(double const *Q) const {
  if (useAdvection) {
    // Coupling is activated!
    const auto Z = NavierStokesSolver_ADERDG_Variables::shortcuts::E + 1;
    return Q[Z];
  }

  // Otherwise, mass fraction of reactant is zero.
  return 0.0;
}

void NavierStokes::PDE::setZ(double *Q, double value) const {
  if (useAdvection) {
    // Coupling is activated!
    const auto Z = NavierStokesSolver_ADERDG_Variables::shortcuts::E + 1;
    Q[Z] = value;
  }
  // Otherwise, nothing happens
}

double NavierStokes::PDE::getHeight(double const *Q) const {
  if (useGravity) {
    // Gravity is included in pressure!
    const auto heightIdx = AbstractNavierStokesSolver_ADERDG::VariableMetrics::SizeVariables;
    return Q[heightIdx];
  }

  // Otherwise, height is zero (-> zero gravitational potential).
  return 0.0;
}

void NavierStokes::PDE::setHeight(double *Q, double value) const {
  if (useGravity) {
    // Gravity is included in pressure!
    const auto heightIdx = AbstractNavierStokesSolver_ADERDG::VariableMetrics::SizeVariables;
    Q[heightIdx] = value;
  }
  // Otherwise, nothing happens
}

void NavierStokes::PDE::setBackgroundState(double *Q, double backgroundRho, double backgroundPressure) const {
  if (useBackgroundState) {
    const auto backgroundStateIdx = AbstractNavierStokesSolver_ADERDG::VariableMetrics::SizeVariables + 1;
      Q[backgroundStateIdx] = backgroundRho;
      Q[backgroundStateIdx+1] = backgroundPressure;
  }
}

std::pair<double, double> NavierStokes::PDE::getBackgroundState(double const *Q) const {
 if (useBackgroundState) {
    const auto backgroundStateIdx = AbstractNavierStokesSolver_ADERDG::VariableMetrics::SizeVariables + 1;
    return {Q[backgroundStateIdx], Q[backgroundStateIdx+1]};
  }
 return {0.0, 0.0};
}

double NavierStokes::PDE::evaluateEnergy(double rho, double pressure, const tarch::la::Vector<DIMENSIONS,double> &j,
        double Z, double height) const {
  const auto invRho = 1./rho;
  const auto chemicalEnergy = q0 * Z;

  // TODO(Lukas) Refactor gravity!
  auto gravityPotential = 0.0;
  if (useGravity) {
    gravityPotential = rho * gravitation * height;
  }

  return pressure/(gamma - 1) + 0.5 * (invRho * j * j) + chemicalEnergy + gravityPotential;
}

double NavierStokes::PDE::evaluatePressure(double E, double rho, const tarch::la::Vector<DIMENSIONS,double> &j,
        double Z, double height) const {
  const auto chemicalEnergy = q0 * Z;

  // TODO(Lukas) Refactor gravity!
  auto gravityPotential = 0.0;
  if (useGravity) {
    gravityPotential = rho * gravitation * height;
  }

  return (gamma-1) * (E - 0.5 * (1.0/rho) * j * j - chemicalEnergy - gravityPotential);
}

double NavierStokes::PDE::evaluateTemperature(double rho, double pressure) const {
  return pressure/(gasConstant * rho);
}

double NavierStokes::PDE::evaluateHeatConductionCoeff(double viscosity) const {
  return 1./Pr * viscosity * gamma * c_v;
}


double NavierStokes::PDE::evaluateViscosity(double T) const {
  return referenceViscosity;
}

void NavierStokes::PDE::evaluateEigenvalues(const double* const Q, const int d, double* lambda) const {
  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j(), getZ(Q), getHeight(Q));
  //assertion5(std::isfinite(p), p, vars.E(), vars.rho(), vars.j(), getZ(Q));

  const double u_n = vars.j(d)/vars.rho();
  const double temperature = evaluateTemperature(vars.rho(), p);
  //const double c = std::sqrt(gamma * gasConstant * temperature);
  const double c = std::sqrt(gamma * (p/vars.rho()));

  //assertion3(std::isfinite(u_n), u_n, vars.j(d), vars.rho());
  //assertion6(std::isfinite(temperature) && temperature >= 0.0, temperature, vars.rho(), vars.E(), vars.j() * vars.j(), p, u_n);
  //assertion3(std::isfinite(c), c, u_n, temperature);

  std::fill_n(lambda, vars.variables(), u_n);
  lambda[0] -= c;
  lambda[1] += c;
}

void NavierStokes::PDE::evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const {
  ReadOnlyVariables vars(Q);

  //const auto pressure = evaluatePressure(vars.E(), vars.rho(), vars.j(), getZ(Q), getHeight(Q));
  //const double pressure = evaluatePressure(vars.E(), vars.rho(), vars.j());

  //const double T = evaluateTemperature(vars.rho(), pressure);
  const double viscosity = referenceViscosity; //evaluateViscosity(T);

  // TODO(Lukas): Need to init to zero here?
  std::fill_n(lambda, vars.variables(), 0.0);
  lambda[0] = (4./3.) * (viscosity/vars.rho());
  lambda[1] = (gamma * viscosity) / (Pr * vars.rho());
  lambda[2] = molecularDiffusionCoeff; // 0 if unused
}

void NavierStokes::PDE::evaluateFlux(const double* Q, const double* gradQ, double** F, bool useViscosity,
				     bool reconstructGradT, double reconstructedGradT,
				     bool gradContainsParameters) const {
  // Variable shortcuts
  const auto rho = NavierStokesSolver_ADERDG_Variables::shortcuts::rho;
  const auto j = NavierStokesSolver_ADERDG_Variables::shortcuts::j;
  const auto E = NavierStokesSolver_ADERDG_Variables::shortcuts::E;
  const auto Z = E + 1; // Only defined if coupling is used!

  ReadOnlyVariables vars(Q);

  auto idxF = kernels::idx2(vars.SizeVariables, DIMENSIONS);

  // gradQ contains params for FV solver but not for DG.
  auto sizeGrad = vars.SizeVariables;
  if (gradContainsParameters) {
    sizeGrad += vars.SizeParameters;
  }
  auto idxGradQ = kernels::idx2(DIMENSIONS, sizeGrad);

  const auto invRho = 1/vars.rho(); // Q(1)/(Q(1)*Q(1)+epsilon)

  //p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
  const auto p = evaluatePressure(vars.E(), vars.rho(), vars.j(), getZ(Q), getHeight(Q));
  assertion5(std::isfinite(p), p, vars.E(), vars.rho(), vars.j(), getHeight(Q));

  auto backgroundRho = 0.0;
  auto backgroundPressure = 0.0;
  std::tie(backgroundRho, backgroundPressure) = getBackgroundState(Q);

  double* f = F[0];
  double* g = F[1];
#if DIMENSIONS == 3
  double* h = F[2];
#endif

  /*
  f(1) = Q(2) !- EQN%mu*gradQ(1,1)
  f(2) = irho*Q(2)*Q(2) + p
  f(3) = irho*Q(2)*Q(3)
  f(4) = irho*Q(2)*Q(4)
  f(5) = irho*Q(2)*(Q(5)+p)
   */

  f[rho] = Q[j];
  f[j] = invRho * Q[j] * Q[j] + (p - backgroundPressure);
  f[j+1] = invRho * Q[j] * Q[j+1];
#if DIMENSIONS == 3
  f[j+2] = invRho * Q[j] * Q[j+2];
#endif
  f[E] = invRho * Q[j] * (Q[E] + p);

  /*
  g(1) = Q(3) !- EQN%mu*gradQ(1,2)
  g(2) = irho*Q(3)*Q(2)
  g(3) = irho*Q(3)*Q(3) + p
  g(4) = irho*Q(3)*Q(4)
  g(5) = irho*Q(3)*(Q(5)+p)
   */

  g[rho] = Q[j+1];
  g[j] = invRho * Q[j+1] * Q[j];
  g[j+1] = invRho * Q[j+1] * Q[j+1] + (p - backgroundPressure);
#if DIMENSIONS == 3
  g[j+2] = invRho * Q[j+1] * Q[j+2];
#endif
  g[E] = invRho * Q[j+1] * (Q[E]+p);

  /*
  h(1) = Q(4) !- EQN%mu*gradQ(1,3)
  h(2) = irho*Q(4)*Q(2)
  h(3) = irho*Q(4)*Q(3)
  h(4) = irho*Q(4)*Q(4) + p
  h(5) = irho*Q(4)*(Q(5)+p)
  */

#if DIMENSIONS == 3
  // TODO
  h[rho] = Q[j+2];
  h[j] = invRho * Q[j+2] * Q[j];
  h[j+1] = invRho * Q[j+2] * Q[j+1];
  h[j+2] = invRho * Q[j+2] * Q[j+2] + (p - backgroundPressure);
  h[E] = invRho * Q[j+2] * (Q[E]+p);
#endif


  // Advection-Diffusion
  if (useAdvection) {
    f[Z] = Q[j+0] * Q[Z] / Q[rho];
    g[Z] = Q[j+1] * Q[Z] / Q[rho];
#if DIMENSIONS == 3
    h[Z] = Q[j+2] * Q[Z] / Q[rho];
#endif
  }

     // Viscous:
  /*
 iRho  = 1./Q(1)
 uu    = Q(2)*iRho
 vv    = Q(3)*iRho
 ww    = Q(4)*iRho
   */
  const auto uu = invRho * Q[j];
  const auto vv = invRho * Q[j+1];
#if DIMENSIONS == 3
  const auto ww = invRho * Q[j+2];
#endif

  /*
  mu    = EQN%mu
  kappa = EQN%kappa
  */
  // TODO(Lukas) Support different visc. models?
  const auto mu = referenceViscosity;
  const auto kappa = evaluateHeatConductionCoeff(mu);

  /*
  uux  = iRho*( gradQ(2,1) - uu*gradQ(1,1) )
  vvx  = iRho*( gradQ(3,1) - vv*gradQ(1,1) )
  wwx  = iRho*( gradQ(4,1) - ww*gradQ(1,1) )
  uuy  = iRho*( gradQ(2,2) - uu*gradQ(1,2) )
  vvy  = iRho*( gradQ(3,2) - vv*gradQ(1,2) )
  wwy  = iRho*( gradQ(4,2) - ww*gradQ(1,2) )
  uuz  = iRho*( gradQ(2,3) - uu*gradQ(1,3) )
  vvz  = iRho*( gradQ(3,3) - vv*gradQ(1,3) )
  wwz  = iRho*( gradQ(4,3) - ww*gradQ(1,3) )
  */

  // Derivatives of velocities:
  const auto uux  = invRho * (gradQ[idxGradQ(0,j+0)] - uu * gradQ[idxGradQ(0, rho)]);
  const auto vvx  = invRho * (gradQ[idxGradQ(0,j+1)] - vv * gradQ[idxGradQ(0, rho)]);
#if DIMENSIONS == 3
  const auto wwx  = invRho * (gradQ[idxGradQ(0,j+2)] - ww * gradQ[idxGradQ(0, rho)]);
#endif

  const auto uuy  = invRho * (gradQ[idxGradQ(1,j+0)] - uu * gradQ[idxGradQ(1, rho)]);
  const auto vvy  = invRho * (gradQ[idxGradQ(1,j+1)] - vv * gradQ[idxGradQ(1, rho)]);
#if DIMENSIONS == 3
  const auto wwy  = invRho * (gradQ[idxGradQ(1,j+2)] - ww * gradQ[idxGradQ(1, rho)]);
#endif

#if DIMENSIONS == 3
  const auto uuz  = invRho * (gradQ[idxGradQ(2,j+0)] - uu * gradQ[idxGradQ(2, rho)]);
  const auto vvz  = invRho * (gradQ[idxGradQ(2,j+1)] - vv * gradQ[idxGradQ(2, rho)]);
  const auto wwz  = invRho * (gradQ[idxGradQ(2,j+2)] - ww * gradQ[idxGradQ(2, rho)]);
#endif

  /*
  icv   = 1./EQN%cv
  iRho2 = iRho*iRho
  iRho3 = iRho2*iRho
  */
  const auto invCv = 1/c_v;
  const auto invRho2 = invRho * invRho;
  const auto invRho3 = invRho2 * invRho;

  /*
  dTdW1 = - Q(5)*iRho2 + iRho3*( Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4) )
  dTdW2 = - Q(2)*iRho2
  dTdW3 = - Q(3)*iRho2
  dTdW4 = - Q(4)*iRho2
  */
#if DIMENSIONS == 2
  const auto dTdW1 = -1 * Q[E] * invRho2 +
          invRho3 * (Q[j]*Q[j] + Q[j+1] * Q[j+1]);
#else
  const auto dTdW1 = -1 * Q[E] * invRho2 +
          invRho3 * (Q[j]*Q[j] + Q[j+1] * Q[j+1] + Q[j+2] * Q[j+2]);
#endif
  const auto dTdW2 = -1 * Q[j] * invRho2;
  const auto dTdW3 = -1 * Q[j+1] * invRho2;
#if DIMENSIONS == 3
  const auto dTdW4 = -1 * Q[j+2] * invRho2;
#endif

   auto Tx = 0.0;
   auto Ty = 0.0;
   auto Tz = 0.0;

   if (reconstructGradT) {
     // Compute derivative from hydrostatic equilibrium.
#if DIMENSIONS == 2
     Ty = reconstructedGradT;
#else
     Tz = reconstructedGradT;
#endif
   } else {
     // Compute derivative from conserved variables.

     /*
     Tx = icv*( dTdW1*gradQ(1,1) + dTdW2*gradQ(2,1) + dTdW3*gradQ(3,1) + dTdW4*gradQ(4,1) + iRho*gradQ(5,1) )
     Ty = icv*( dTdW1*gradQ(1,2) + dTdW2*gradQ(2,2) + dTdW3*gradQ(3,2) + dTdW4*gradQ(4,2) + iRho*gradQ(5,2) )
     Tz = icv*( dTdW1*gradQ(1,3) + dTdW2*gradQ(2,3) + dTdW3*gradQ(3,3) + dTdW4*gradQ(4,3) + iRho*gradQ(5,3) )
     */
#if DIMENSIONS == 2
     Tx = invCv * (dTdW1 * gradQ[idxGradQ(0, rho)] +
          dTdW2 * gradQ[idxGradQ(0, j)] +
          dTdW3 * gradQ[idxGradQ(0, j+1)] +
          invRho * gradQ[idxGradQ(0,E)]);

     Ty = invCv * (dTdW1 * gradQ[idxGradQ(1, rho)] +
          dTdW2 * gradQ[idxGradQ(1, j)] +
          dTdW3 * gradQ[idxGradQ(1, j+1)] +
          invRho * gradQ[idxGradQ(1,E)]);
#else
    // TODO(Lukas): Untested for 3D!
     Tx = invCv * (dTdW1 * gradQ[idxGradQ(0, rho)] +
          dTdW2 * gradQ[idxGradQ(0, j)] +
          dTdW3 * gradQ[idxGradQ(0, j+1)] +
          dTdW4 * gradQ[idxGradQ(0, j+2)] +
          invRho * gradQ[idxGradQ(0,E)]);

     Ty = invCv * (dTdW1 * gradQ[idxGradQ(1, rho)] +
          dTdW2 * gradQ[idxGradQ(1, j)] +
          dTdW3 * gradQ[idxGradQ(1, j+1)] +
          dTdW4 * gradQ[idxGradQ(1, j+2)] +
          invRho * gradQ[idxGradQ(1,E)]);

     Tz = invCv * (dTdW1 * gradQ[idxGradQ(2, rho)] +
          dTdW2 * gradQ[idxGradQ(2, j)] +
          dTdW3 * gradQ[idxGradQ(2, j+1)] +
          dTdW4 * gradQ[idxGradQ(2, j+2)] +
          invRho * gradQ[idxGradQ(2,E)]);
#endif
   }
  /*
  divV23  = 2./3.*(uux + vvy + wwz)
  */
#if DIMENSIONS == 2
  const auto divV23 = 2./3. * (uux + vvy);
#else
  const auto divV23 = 2./3. * (uux + vvy + wwz);
#endif

  // Arrays of viscous fluxes.
  // These do NOT contain temperature diffusion terms!
  double Fv[vars.Size] = {0.0};
  double Gv[vars.Size] = {0.0};
  double Hv[vars.Size] = {0.0};
  /*
  Fv(1) = 0.
  Fv(2) = mu*( 2*uux - divV23 )
  Fv(3) = mu*(   uuy + vvx    )
  Fv(4) = mu*(   uuz + wwx    )
  Fv(5) = Fv(2)*uu + Fv(3)*vv + Fv(4)*ww + kappa*Tx
  */
  Fv[rho] = 0.0;
  Fv[j+0] = mu * (2 * uux - divV23);
  Fv[j+1] = mu * (uuy + vvx);
#if DIMENSIONS == 2
  Fv[E] = Fv[j] * uu + Fv[j+1] * vv;
#else
  Fv[j+2] = mu * (uuz + wwx);
  Fv[E] = Fv[j] * uu + Fv[j+1] * vv + Fv[j+2] * ww;
#endif

  /*
  Gv(1) = 0.
  Gv(2) = Fv(3)
  Gv(3) = mu*( 2*vvy - divV23 )
  Gv(4) = mu*(   vvz + wwy    )
  Gv(5) = Gv(2)*uu + Gv(3)*vv + Gv(4)*ww + kappa*Ty
  */
  Gv[rho] = 0.0;
  Gv[j+0] = Fv[j+1];
  Gv[j+1] = mu * (2 * vvy - divV23);
#if DIMENSIONS == 2
  Gv[E] = Gv[j] * uu + Gv[j+1] * vv;
#else
  Gv[j+2] = mu * (vvz + wwy);
  Gv[E] = Gv[j] * uu + Gv[j+1] * vv + Gv[j+2] * ww;
#endif

  // TODO(Lukas) Test viscous flux for 3D!
  /*
  Hv(1) = 0.
  Hv(2) = Fv(4)
  Hv(3) = Gv(4)
  Hv(4) = mu*( 2*wwz - divV23 )
  Hv(5) = Hv(2)*uu + Hv(3)*vv + Hv(4)*ww + kappa*Tz
  */
#if DIMENSIONS == 3
  Hv[rho] = 0.0;
  Hv[j+0] = Fv[j+2];
  Hv[j+1] = Gv[j+2];
  Hv[j+2] = mu * (2 * wwz - divV23);
  Hv[E] = Hv[j] * uu + Hv[j+1] * vv + Hv[j+2] * ww;
#endif

  /*
  f = f - Fv
  g = g - Gv
  h = h - Hv
  */
  if (useViscosity) {
    for (int i = 0; i < vars.Size; ++i) {
      f[i] -= Fv[i];
      g[i] -= Gv[i];
#if DIMENSIONS == 3
      h[i] -= Hv[i];
#endif
    }
  }

  // Heat flux changes for reactive NS-Equation due to change in
  // definition of pressure.
  // TODO(Lukas) Maybe fix for reconstructed heat flux!
  if (useAdvection) {
    const double factor = invCv * q0 * invRho2;
    Tx -= factor *
            (Q[rho] * gradQ[idxGradQ(0,Z)] -
             Q[Z] * gradQ[idxGradQ(0, rho)]);
    Ty -= factor *
            (Q[rho] * gradQ[idxGradQ(1,Z)] -
             Q[Z] * gradQ[idxGradQ(1, rho)]);
#if DIMENSIONS == 3
    Tz -= factor *
             (Q[rho] * gradQ[idxGradQ(2,Z)] -
              Q[Z] * gradQ[idxGradQ(2, rho)]);
#endif
  }

  // Including the gravitational potential in the pressure
  // of course implies that we need to change the derivative of the temperature as well,
  // but only in z direction
  if (useGravity) {
    // Gravity happens in y direction in 2D and in z direction in 3D
#if DIMENSIONS == 2
    Ty -= (gravitation * (gamma - 1)) / gasConstant;
#else
    Tz -= (gravitation * (gamma - 1)) / gasConstant;
#endif
  }

  // Heat flux
  f[E] -= kappa * Tx;
  g[E] -= kappa * Ty;
# if DIMENSIONS == 3
  h[E] -= kappa * Tz;
#endif

  // Molecular diffusion for Advection-Reaction-Diffusion part.
  if (molecularDiffusionCoeff > 0) {
    f[rho] -= molecularDiffusionCoeff * gradQ[idxGradQ(0, rho)];
    g[rho] -= molecularDiffusionCoeff * gradQ[idxGradQ(1, rho)];

    f[Z] -= molecularDiffusionCoeff * gradQ[idxGradQ(0, Z)];
    g[Z] -= molecularDiffusionCoeff * gradQ[idxGradQ(1, Z)];

#if DIMENSIONS == 3
    h[rho] -= molecularDiffusionCoeff * gradQ[idxGradQ(2, rho)];
    h[Z] -= molecularDiffusionCoeff * gradQ[idxGradQ(2, Z)];
#endif
  }

}


double NavierStokes::PDE::evaluatePotentialTemperature(double temperature, double pressure) const {
    return temperature / std::pow((pressure / referencePressure), (gasConstant / c_p));
}
