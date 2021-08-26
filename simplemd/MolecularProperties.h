// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULARPROPERTIES_H_
#define _MOLECULARDYNAMICS_MOLECULARPROPERTIES_H_

namespace simplemd {
  class MolecularProperties;
}


/** describes the properties of a Lennard-Jones fluid.
 *  @author Philipp Neumann
 */
class simplemd::MolecularProperties {

  public:
    MolecularProperties(
      const double& mass,
      const double& epsilon,
      const double& sigma,
      const double& cutoffRadius,
      const double& kB
    ): _mass(mass), _epsilon(epsilon), _sigma(sigma), _cutoffRadius(cutoffRadius), _kB(kB){}
    ~MolecularProperties(){}

    const double& getMass() const {return _mass;}
    const double& getEpsilon() const {return _epsilon;}
    const double& getSigma() const {return _sigma;}
    const double& getCutOffRadius() const { return _cutoffRadius;}
    const double& getKB() const { return _kB;}

  private:
    /** mass of a fluid molecule */
    double _mass;

    /** epsilon parameter of the LJ-description */
    double _epsilon;

    /** sigma parameter of the LJ-description */
    double _sigma;

    /** cutoff radius */
    double _cutoffRadius;

    /** Boltzmann's constant */
    double _kB;
};

#endif

