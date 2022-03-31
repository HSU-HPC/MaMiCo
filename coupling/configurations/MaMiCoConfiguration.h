// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MAMICOCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MAMICOCONFIGURATION_H_

#include "coupling/configurations/BoundaryForceConfiguration.h"
#include "coupling/configurations/MacroscopicCellConfiguration.h"
#include "coupling/configurations/MomentumInsertionConfiguration.h"
#include "coupling/configurations/NoiseReductionConfiguration.h"
#include "coupling/configurations/ParallelTopologyConfiguration.h"
#include "coupling/configurations/ParticleInsertionConfiguration.h"
#include "coupling/configurations/ThermostatConfiguration.h"
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
        _isDefinedTransferStrategy(false), _isDefinedNoiseReduction(false), _isDefinedParallelTopology(false), _isDefinedThermostat(false) {}

  /** Destructor */
  virtual ~MaMiCoConfiguration() {}

  /** parseSubtag
   * 	@param node
   */
  void parseSubtag(tinyxml2::XMLElement *node);

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "mamico"; }

  /** checks if the configuration is valid. This operation usually fails, if
e.g.
         *	1. parseSubtag() hasn't been called, i.e. configuration has not been used,
or
     *  2. parseSubtag() failed due to a wrong file.
         * 	@return _isValid
     */
  bool isValid() const { return _isValid; }

  /**
   * 	@return _macroscopicCellConfiguration
   */
  const coupling::configurations::MacroscopicCellConfiguration<dim> &getMacroscopicCellConfiguration() const { return _macroscopicCellConfiguration; }
  /**
   * 	@return _particleInsertionConfiguration
   */
  const coupling::configurations::ParticleInsertionConfiguration &getParticleInsertionConfiguration() const {
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
  const coupling::configurations::MomentumInsertionConfiguration &getMomentumInsertionConfiguration() const {
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
  const coupling::configurations::BoundaryForceConfiguration<dim> &getBoundaryForceConfiguration() const {
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
  const coupling::configurations::TransferStrategyConfiguration<dim> &getTransferStrategyConfiguration() const {
    if (!_isDefinedTransferStrategy) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Transfer-Strategy not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _transferStrategyConfiguration;
  }

  /**
   * 	@return _noiseReductionConfiguration
   */
  const coupling::configurations::NoiseReductionConfiguration &getNoiseReductionConfiguration() const {
    if (!_isDefinedNoiseReduction) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Noise-Reduction not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _noiseReductionConfiguration;
  }

  /**
   * 	@return _parallelTopologyConfiguration
   */
  const coupling::configurations::ParallelTopologyConfiguration &getParallelTopologyConfiguration() const {
    if (!_isDefinedParallelTopology) {
      std::cout << "ERROR coupling::configurations::MaMiCoConfiguration: "
                   "Parallel-Topology not defined!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    return _parallelTopologyConfiguration;
  }

  const coupling::configurations::ThermostatConfiguration &getThermostatConfiguration() const {
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
  coupling::configurations::MacroscopicCellConfiguration<dim> _macroscopicCellConfiguration;
  coupling::configurations::ParticleInsertionConfiguration _particleInsertionConfiguration;
  coupling::configurations::MomentumInsertionConfiguration _momentumInsertionConfiguration;
  coupling::configurations::BoundaryForceConfiguration<dim> _boundaryForceConfiguration;
  coupling::configurations::TransferStrategyConfiguration<dim> _transferStrategyConfiguration;
  coupling::configurations::NoiseReductionConfiguration _noiseReductionConfiguration;
  coupling::configurations::ParallelTopologyConfiguration _parallelTopologyConfiguration;
  coupling::configurations::ThermostatConfiguration _thermostatConfiguration;
  bool _isDefinedParticleInsertion;
  bool _isDefinedMomentumInsertion;
  bool _isDefinedBoundaryForce;
  bool _isDefinedTransferStrategy;
  bool _isDefinedNoiseReduction;
  bool _isDefinedParallelTopology;
  bool _isDefinedThermostat;
};
#include "coupling/configurations/MaMiCoConfiguration.cpph"

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MACROSCOPICCELLCONFIGURATION_H_
