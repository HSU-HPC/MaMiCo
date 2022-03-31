#ifndef _MOLECULARDYNAMICS_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_
#define _MOLECULARDYNAMICS_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_

#include "tarch/utils/RandomNumberService.h"

#include "simplemd/LinkedCell.h"
#include "simplemd/MolecularDynamicsDefinitions.h"

namespace simplemd {
namespace cellmappings {
class VaryCheckpointMapping;
}
} // namespace simplemd

class simplemd::cellmappings::VaryCheckpointMapping {
public:
  VaryCheckpointMapping(const double &molecularMass, const double &kB, const double &temperature, const double &sigma,
                        const tarch::la::Vector<MD_DIM, double> &meshSize)
      : _molecularMass(molecularMass), _kB(kB), _temperature(temperature), _sigma(sigma), _meshSize(meshSize) {}

  void beginCellIteration() {}

  void endCellIteration() {}

  void handleCell(LinkedCell &cell, const unsigned int &cellIndex) {
    tarch::la::Vector<MD_DIM, double> meanVelocityForCell(0.0);
    double stdDeviation = std::sqrt(MD_DIM * _kB * _temperature / _molecularMass);

    tarch::la::Vector<MD_DIM, double> randomNumbers(0.0);

    for (std::list<Molecule *>::iterator molecule = cell.begin(); molecule != cell.end(); molecule++) {
      meanVelocityForCell += (*molecule)->getVelocity();
    }
    meanVelocityForCell = meanVelocityForCell / (double)cell.getConstList().size() / _molecularMass;

    for (std::list<Molecule *>::iterator molecule = cell.begin(); molecule != cell.end(); molecule++) {
      randomNumbers[0] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      for (unsigned int d = 1; d < MD_DIM; ++d) {
        randomNumbers[d] = tarch::utils::RandomNumberService::getInstance().getGaussianRandomNumber();
      }

      tarch::la::Vector<MD_DIM, double> &mVelocity = (*molecule)->getVelocity();
#if (MD_DIM == 1)
      mVelocity = meanVelocityForCell + stdDeviation * randomNumbers;
#elif (MD_DIM == 2)
      mVelocity[0] = meanVelocityForCell[0] + stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
      mVelocity[1] = meanVelocityForCell[1] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]));
#elif (MD_DIM == 3)
      mVelocity[0] = meanVelocityForCell[0] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::cos(randomNumbers[2]));
      mVelocity[1] = meanVelocityForCell[1] + stdDeviation * (randomNumbers[0] * std::sin(randomNumbers[1]) * std::sin(randomNumbers[2]));
      mVelocity[2] = meanVelocityForCell[2] + stdDeviation * (randomNumbers[0] * std::cos(randomNumbers[1]));
#endif

      tarch::la::Vector<MD_DIM, double> &mPosition = (*molecule)->getPosition();
      mPosition = mPosition + 1e-6 * mVelocity;
    }
  }

private:
  const double _molecularMass;
  const double _kB;
  const double _temperature;
  const double _sigma;
  const tarch::la::Vector<MD_DIM, double> _meshSize;
};

#endif //_MOLECULARDYNAMICS_CELLMAPPINGS_VARYCHECKPOINTMAPPING_H_