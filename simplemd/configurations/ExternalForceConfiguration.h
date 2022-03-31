// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_EXTERNALFORCECONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_EXTERNALFORCECONFIGURATION_H_

#include "simplemd/MolecularDynamicsUserInput.h"
#include "tarch/la/Vector.h"
#include "tarch/configuration/Configuration.h"

namespace simplemd {
namespace configurations { class ExternalForceConfiguration; }
}

/** parses a constant, external force term.
 *  @author Philipp Neumann
 */
class simplemd::configurations::ExternalForceConfiguration
    : public tarch::configuration::Configuration {
public:
  ExternalForceConfiguration() : _externalForce(0.0), _isValid(false) {}
  virtual ~ExternalForceConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement *node);

  /**
   * Return name of xml tag that is associated to the configuration.
   */
  std::string getTag() const { return "external-force"; }

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
  bool isValid() const { return _isValid; }

  /** getters for all parsed and computed quantities */
  const tarch::la::Vector<MD_DIM, double> &getExternalForce() const {
    return _externalForce;
  }

private:
  static const std::string VALUE;

  tarch::la::Vector<MD_DIM, double> _externalForce;
  bool _isValid;
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_EXTERNALFORCECONFIGURATION_H_
