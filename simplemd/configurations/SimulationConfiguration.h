// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_SIMULATIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_SIMULATIONCONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace simplemd {
namespace configurations {
class SimulationConfiguration;
}
} // namespace simplemd

/** simulation configuration for MD. Currently includes information on timestepping only.
 *  @author Philipp Neumann
 */
class simplemd::configurations::SimulationConfiguration : public tarch::configuration::Configuration {
public:
  SimulationConfiguration();
  virtual ~SimulationConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement* node);

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
  const double& getDt() const { return _dt; }
  const unsigned int& getNumberOfTimesteps() const { return _numberOfTimesteps; }
  const unsigned int& computeMacroscopicQuantitiesEveryTimestep() const { return _computeMacroscopicQuantitiesEveryTimestep; }
  const bool& fixSeed() const { return _fixSeed; }
  const bool& useOverlappingCommunicationWithForceComputation() const { return _overlapCommWithForceComputation; }

private:
  static const std::string DT;
  static const std::string NUMBER_OF_TIMESTEPS;
  static const std::string COMPUTE_QUANTITIES_EVERY_TIMESTEP;
  static const std::string FIX_SEED;
  static const std::string OVERLAP_COMMUNICATION_WITH_FORCE_COMPUTATION;

  double _dt;
  unsigned int _numberOfTimesteps;

  /** number of timesteps between subsequent macroscopic quantity evaluations */
  unsigned int _computeMacroscopicQuantitiesEveryTimestep;

  /** if true, the seed of the random number service is fixed in all simulation runs. */
  bool _fixSeed;

  /** if true, the force computation and the send/receive-operations for molecules at the process boundaries
   *  are overlapped.
   */
  bool _overlapCommWithForceComputation;

  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_SIMULATIONCONFIGURATION_H_
