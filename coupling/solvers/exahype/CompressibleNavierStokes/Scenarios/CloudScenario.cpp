#include "CloudScenario.h"
#include "Atmosphere.h"

NavierStokes::CloudScenario::Bubble::Bubble(
    BubbleType bubbleType, const double tempDifference, const double size,
    const double decay, const double centerX, const double centerZ)
    : bubbleType(bubbleType),
      tempDifference(tempDifference),
      size(size),
      decay(decay),
      centerX(centerX),
      centerY(centerX), // TODO(Lukas) extend this.
      centerZ(centerZ) {}

double NavierStokes::CloudScenario::Bubble::evaluatePerturbation(
			    double posX, double posY, double posZ) const {
  auto perturbation = 0.0;
  const auto distX = posX - centerX;
#if DIMENSIONS == 2
  const auto distY = 0.0;
#else
  const auto distY = posY - centerY;
#endif
  const auto distZ = posZ - centerZ;
  const auto distanceToCenter = std::sqrt(distX * distX +  distY * distY + distZ * distZ);

  if (bubbleType == BubbleType::smooth) {
    if (distanceToCenter <= size) {
      perturbation = tempDifference;

    } else {
      const auto d = distanceToCenter - size;
      perturbation = tempDifference * std::exp(-(d * d) / (decay * decay));
    }
  } else {
    // bubbleType == BubbleType::cosine
    if (distanceToCenter <= size) {
      const auto pi = std::acos(-1);  // not constexpr for Intel compiler!
      perturbation =
          tempDifference / 2 * (1 + std::cos((pi * distanceToCenter) / size));
    }
  }
  assertion3(std::isfinite(perturbation), posX, posY, posZ);
  return perturbation;
}

void NavierStokes::CloudScenario::initialValues(const double* const x,
                                                const PDE& ns,
                                                Variables& vars) {
  initialValues(x, ns, vars, 0.0);
}

void NavierStokes::CloudScenario::initialValues(const double* const x,
                                                const PDE& ns, Variables& vars,
                                                double initialZ) {
  const auto posX = x[0];
  const auto posY = (DIMENSIONS == 2) ? 0.0 : x[1];
  const auto posZ = (DIMENSIONS == 2) ? x[1] : x[2];

  double potentialT = getBackgroundPotentialTemperature();

  double Z = 0.0;  // Only used for coupling test!

  for (const auto& bubble : bubbles) {
    potentialT += bubble.evaluatePerturbation(posX, posY, posZ);
    Z += bubble.evaluatePerturbation(posX, posY, posZ);  // TODO(Lukas): Fix value.
  }

  // Air is initially at rest.
#if DIMENSIONS == 2
  vars.j(0, 0);
#elif DIMENSIONS == 3
  vars.j(0, 0, 0);
#endif

  // First compute overall state
  const auto pressure = computeHydrostaticPressure(
      ns, getGravity(), posZ, getBackgroundPotentialTemperature());
  const auto temperature = potentialTToT(ns, pressure, potentialT);
  vars.rho() = pressure / (ns.gasConstant * temperature);
  vars.E() = ns.evaluateEnergy(vars.rho(), pressure, vars.j(), Z,
                               ns.getHeight(vars.data()));
  assertion4(std::isfinite(vars.rho()), vars.rho(), vars.E(), temperature,
             pressure);
  assertion4(std::isfinite(vars.E()), vars.rho(), vars.E(), temperature,
             pressure);
  ns.setZ(vars.data(), Z);

  if (ns.useBackgroundState) {
    // Then compute background state (without pot.T. perturbation
    const auto backgroundTemperature =
        potentialTToT(ns, pressure, getBackgroundPotentialTemperature());
    const auto backgroundRho =
        pressure / (ns.gasConstant * backgroundTemperature);
    ns.setBackgroundState(vars.data(), backgroundRho, pressure);
  }
}

double NavierStokes::CloudScenario::getGamma() const { return gamma; }

double NavierStokes::CloudScenario::getPr() const { return Pr; }

double NavierStokes::CloudScenario::getC_v() const { return c_v; }

double NavierStokes::CloudScenario::getC_p() const { return c_p; }

double NavierStokes::CloudScenario::getGasConstant() const {
  return gasConstant;
}

double NavierStokes::CloudScenario::getReferencePressure() const {
  return referencePressure;
}

double NavierStokes::CloudScenario::getGravity() const { return 9.81; }

double NavierStokes::CloudScenario::getBackgroundPotentialTemperature() const {
  return 300;
}

NavierStokes::BoundaryType NavierStokes::CloudScenario::getBoundaryType(
    int faceId) {
  return BoundaryType::hydrostaticWall;
}

void NavierStokes::CloudScenario::source(
    const tarch::la::Vector<DIMENSIONS, double>& x, double t, const PDE& ns,
    const double* const Q, double* S) {
  Scenario::source(x, t, ns, Q, S);

  double rhoPertubation = Q[0];
  if (ns.useBackgroundState) {
    auto backgroundRho = 0.0;
    auto backgroundPressure = 0.0;
    std::tie(backgroundRho, backgroundPressure) = ns.getBackgroundState(Q);

    rhoPertubation -= backgroundRho;
  }
  S[DIMENSIONS] = -1 * rhoPertubation * getGravity();

  // Only use this source term if the gravitational force is not already
  // included in the pressure.
  if (!ns.useGravity) {
    S[DIMENSIONS + 1] = -1 * Q[2] * getGravity();
  }
}
