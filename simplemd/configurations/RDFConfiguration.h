// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_RDFCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_RDFCONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace simplemd {
namespace configurations {
class RDFConfiguration;
}
} // namespace simplemd

/** configuration input for RDF sampling.
 *  @author Philipp Neumann
 */
class simplemd::configurations::RDFConfiguration : public tarch::configuration::Configuration {
public:
  RDFConfiguration();
  virtual ~RDFConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement *node);

  /**
   * Return name of xml tag that is associated to the configuration.
   */
  virtual std::string getTag() const;

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
  const bool &isDefined() const { return _isDefined; }
  const unsigned int &getStartAtTimestep() const { return _startAtTimestep; }
  const unsigned int &getEvaluateEveryTimestep() const { return _evaluateEveryTimestep; }
  const unsigned int &getWriteEveryTimestep() const { return _writeEveryTimestep; }
  const unsigned int &getNumberOfPoints() const { return _numberPoints; }

private:
  static const std::string START_AT_TIMESTEP;
  static const std::string EVALUATE_EVERY_TIMESTEP;
  static const std::string WRITE_EVERY_TIMESTEP;
  static const std::string NUMBER_OF_POINTS;

  /** first timestep when the rdf-evaluation shall be started */
  unsigned int _startAtTimestep;

  /** timesteps per subsequent rdf evaluations */
  unsigned int _evaluateEveryTimestep;

  /** number of timesteps per vtk output. If this value is zero, no output is
   * written at all */
  unsigned int _writeEveryTimestep;

  /** number of points used to discretise the distance between two molecules
   * (between r=0 and r=r_cutoff) */
  unsigned int _numberPoints;

  /** true, if the configuration was parsed */
  bool _isDefined;

  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_RDFCONFIGURATION_H_
