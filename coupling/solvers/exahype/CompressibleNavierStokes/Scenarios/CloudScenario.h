#ifndef COMPRESSIBLENAVIERSTOKES_CLOUDSCENARIO_H
#define COMPRESSIBLENAVIERSTOKES_CLOUDSCENARIO_H

#include "Scenario.h"
namespace NavierStokes {

class CloudScenario : public Scenario {
 public:
  enum class BubbleType { smooth, cosine };

  struct Bubble {
    Bubble(BubbleType bubbleType, const double tempDifference,
           const double size, const double decay, const double centerX,
           const double centerZ);

    // temp. difference
    const double tempDifference;
    const double size;     // [m]
    const double decay;    // [m]
    const double centerX;  // [m]
    const double centerY;  // [m], currently always equal to centerX
    const double centerZ;  // [m]
    const BubbleType bubbleType;

    double evaluatePerturbation(double posX, double posY, double posZ) const;
  };

  void initialValues(const double* const x, const PDE& ns,
                     Variables& vars) override;
  void initialValues(const double* const x, const PDE& ns, Variables& vars,
                     double initialZ);

  void source(const tarch::la::Vector<DIMENSIONS, double>& x, double t,
              const PDE& ns, const double* const Q, double* S) override;

  BoundaryType getBoundaryType(int faceId);

  std::vector<Bubble> bubbles;

  // Constants for dry air.
  const double gamma = 1.4;
  const double Pr = 0.71;
  const double gasConstant = 287.058;
  // const double c_p = 1.005 * 1000;
  const double c_p = gamma / (gamma - 1) * gasConstant;
  const double c_v = 1 / (gamma - 1) * gasConstant;
  const double referencePressure = 10000;

  double getGamma() const override;
  double getPr() const override;
  double getC_v() const override;
  double getC_p() const override;
  double getGasConstant() const override;
  double getReferencePressure() const override;
  double getGravity() const override;
  double getBackgroundPotentialTemperature() const override;
};

}  // namespace NavierStokes

#endif  // COMPRESSIBLENAVIERSTOKES_CLOUDSCENARIO_H
