#ifndef NAVIERSTOKES_PDE_H
#define NAVIERSTOKES_PDE_H

#include <cmath>
#include <memory>
#include "tarch/la/Vector.h"
#include "kernels/KernelUtils.h"
#include "NavierStokesSolver_ADERDG_Variables.h"

namespace NavierStokes {
  class PDE;
  class Scenario;
}

using Variables = NavierStokes::AbstractNavierStokesSolver_ADERDG::Variables;
using ReadOnlyVariables = NavierStokes::AbstractNavierStokesSolver_ADERDG::ReadOnlyVariables;
using Fluxes = NavierStokes::AbstractNavierStokesSolver_ADERDG::Fluxes;

class NavierStokes::PDE {
public:
  PDE();
  PDE(double referenceViscosity, double referencePressure, double gamma, double Pr, double c_v,
          double c_p, double gasConstant);
  PDE(double referenceViscosity, Scenario& scenario, bool useGravity, bool useBackgroundState);

  // Getter/Setters for optional variables/parameters.
  double getZ(double const *Q) const;
  void setZ(double *Q, double value) const;
  double getHeight(double const *Q) const;
  void setHeight(double *Q, double value) const;
  void setBackgroundState(double *Q, double backgroundRho, double backgroundPressure) const;
  std::pair<double, double> getBackgroundState(double const *Q) const;

  double evaluateEnergy(double rho, double pressure, const tarch::la::Vector<DIMENSIONS,double> &j,
          double Z=0.0, double height=0.0) const;
  double evaluatePressure(double E, double rho, const tarch::la::Vector<DIMENSIONS,double> &j,
          double Z=0.0, double height=0.0) const;
  double evaluateTemperature(double rho, double pressure) const;
  double evaluatePotentialTemperature(double temperature, double pressure) const;
  double evaluateHeatConductionCoeff(double viscosity) const;  
  double evaluateViscosity(double T) const;
  void evaluateEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateDiffusiveEigenvalues(const double* const Q, const int d, double* lambda) const;
  void evaluateFlux(const double* Q, const double* gradQ, double** F, bool useViscosity=true,
		    bool reconstructGradT=false, double reconstructedGradT=0.0,
		    bool gradContainsParameters=false) const;

  double referenceViscosity; // [Pa s] = [kg/(ms)] Reference dynamic viscosity
  double referencePressure; // [Pa], reference pressure used by some scenarios.

  double gamma; // [1] Ratio of specific heats, c_p/c_v
  double Pr; // [1] Prandtl number

  // [m^2/(s^2 K)] = [J/(kg K)]
  double c_v; // [J/(kg K)] Specific heat capacity at constant volume
  double c_p; // [J/(kg K)] Specific heat capacity at constant pressure
  double gasConstant; // [J/(kg K)] Specific gas constant

  // Advection-Diffusion Part
  double q0; // ???
  double molecularDiffusionCoeff; // ???
  bool useAdvection;

  // Settings for scenarios with gravity source term
  bool useGravity;
  bool useBackgroundState;
  double gravitation;

  // Unused:
  double referenceT;
  double sutherlandC;
  double sutherlandLambda;
};

#endif // NAVIERSTOKES_PDE_H
