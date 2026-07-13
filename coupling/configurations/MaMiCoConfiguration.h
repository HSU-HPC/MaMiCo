// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MAMICOCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MAMICOCONFIGURATION_H_

#include "coupling/configurations/BoundaryForceConfiguration.h"
#include "coupling/configurations/CouplingCellConfiguration.h"
#include "coupling/configurations/MomentumInsertionConfiguration.h"
#include "coupling/configurations/ParallelTopologyConfiguration.h"
#include "coupling/configurations/ParticleInsertionConfiguration.h"
#include "coupling/configurations/ThermostatConfiguration.h"
#include "coupling/configurations/TimeIntegrationConfiguration.h"
#include "coupling/configurations/TransferStrategyConfiguration.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace coupling {
namespace configurations {
template <unsigned int dim> class MaMiCoConfiguration;
}
} // namespace coupling

/** parses all sub-tags for MaMiCo configuration. Derive from the class
 * tarch::configuration::Configuration
 *	@brief parses all sub-tags for MaMiCo configuration.
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::configurations::MaMiCoConfiguration : public tarch::configuration::Configuration {
public:
  /** Constructor, initializes the class  */
  MaMiCoConfiguration()
      : _isValid(true), _isDefinedParticleInsertion(false), _isDefinedMomentumInsertion(false), _isDefinedBoundaryForce(false),
        _isDefinedTransferStrategy(false), _isDefinedParallelTopology(false), _isDefinedThermostat(false) {}

  /** Destructor */
  virtual ~MaMiCoConfiguration() {}

  /** parseSubtag
   *  @param node
   */
  void parseSubtag(tinyxml2::XMLElement* node);

  /** Returns name of xml tag that is associated to the configuration.
   *  @return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "mamico"; }

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
   *    @return _isValid
   */
  bool isValid() const { return _isValid; }

  /**
   *  @return _couplingCellConfiguration
   */
  const coupling::configurations::CouplingCellConfiguration<dim>& getCouplingCellConfiguration() const { return _couplingCellConfiguration; }

  /**
   *  @return _particleInsertionConfiguration
   */
  const coupling::configurations::ParticleInsertionConfiguration& getParticleInsertionConfiguration() const {
    if (!_isDefinedParticleInsertion) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Particle insertion not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _particleInsertionConfiguration;
  }

  /**
   * 	@return _momentumInsertionConfiguration
   */
  const coupling::configurations::MomentumInsertionConfiguration<dim>& getMomentumInsertionConfiguration() const {
    if (!_isDefinedMomentumInsertion) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Momentum insertion not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _momentumInsertionConfiguration;
  }

  /**
   * 	@return _boundaryForceConfiguration
   */
  const coupling::configurations::BoundaryForceConfiguration<dim>& getBoundaryForceConfiguration() const {
    if (!_isDefinedBoundaryForce) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Boundary force not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _boundaryForceConfiguration;
  }

  /**
   * 	@return _transferStrategyConfiguration
   */
  const coupling::configurations::TransferStrategyConfiguration<dim>& getTransferStrategyConfiguration() const {
    if (!_isDefinedTransferStrategy) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Transfer-Strategy not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _transferStrategyConfiguration;
  }

  /**
   *  @return _parallelTopologyConfiguration
   */
  const coupling::configurations::ParallelTopologyConfiguration& getParallelTopologyConfiguration() const {
    if (!_isDefinedParallelTopology) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Parallel-Topology not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _parallelTopologyConfiguration;
  }

  /**
   *  @return _timeIntegrationConfiguration
   */
  const coupling::configurations::TimeIntegrationConfiguration& getTimeIntegrationConfiguration() const {
    // is optional, thus always defined
    return _timeIntegrationConfiguration;
  }

  /**
   * @return _thermostatConfiguration
   */
  const coupling::configurations::ThermostatConfiguration& getThermostatConfiguration() const {
    if (!_isDefinedThermostat) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Thermostat not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _thermostatConfiguration;
  }

private:
  bool _isValid;
  coupling::configurations::CouplingCellConfiguration<dim> _couplingCellConfiguration;
  coupling::configurations::ParticleInsertionConfiguration _particleInsertionConfiguration;
  coupling::configurations::MomentumInsertionConfiguration<dim> _momentumInsertionConfiguration;
  coupling::configurations::BoundaryForceConfiguration<dim> _boundaryForceConfiguration;
  coupling::configurations::TransferStrategyConfiguration<dim> _transferStrategyConfiguration;
  coupling::configurations::ParallelTopologyConfiguration _parallelTopologyConfiguration;
  coupling::configurations::TimeIntegrationConfiguration _timeIntegrationConfiguration;
  coupling::configurations::ThermostatConfiguration _thermostatConfiguration;
  bool _isDefinedParticleInsertion;
  bool _isDefinedMomentumInsertion;
  bool _isDefinedBoundaryForce;
  bool _isDefinedTransferStrategy;
  bool _isDefinedParallelTopology;
  bool _isDefinedThermostat;
};
#include "coupling/configurations/MaMiCoConfiguration.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MAMICOCONFIGURATION_H_
