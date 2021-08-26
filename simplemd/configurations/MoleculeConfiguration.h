// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULECONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULECONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace simplemd {
  namespace configurations {
    class MoleculeConfiguration;
  }
}


/** reads properties for a single molecule.
 *  @author Philipp Neumann
 */
class simplemd::configurations::MoleculeConfiguration: public tarch::configuration::Configuration {
  public:
    MoleculeConfiguration();
    virtual ~MoleculeConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node );

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const;

    /**
     * Is config valid?
     *
     * This operation usually fails, if
     *
     * - parseSubtag() hasn't been called, i.e. configuration has not been
     *   used, or
     * - parseSubtag() failed due to a wrong file.
     *
     * If a tag ain't optional and parseSubtag() was not called (first case)
     */
    bool isValid() const;

    /** getters for all parsed and computed quantities */
    const tarch::la::Vector<MD_DIM,double>& getMeanVelocity() const { return _meanVelocity;}
    const double& getTemperature() const { return _temperature; }
    const double& getMass() const { return _mass;}
    const double& getEpsilon() const { return _epsilon; }
    const double& getSigma() const { return _sigma; }

  private:
    static const std::string MEAN_VELOCITY;
    static const std::string TEMPERATURE;
    static const std::string MASS;
    static const std::string EPSILON;
    static const std::string SIGMA;

    tarch::la::Vector<MD_DIM,double> _meanVelocity;
    double _temperature;
    double _mass;
    double _epsilon;
    double _sigma;

    bool _isValid;
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULECONFIGURATION_H_
