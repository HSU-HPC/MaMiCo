// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULARDYNAMICSCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULARDYNAMICSCONFIGURATION_H_

#include "simplemd/configurations/CheckpointConfiguration.h"
#include "simplemd/configurations/DomainConfiguration.h"
#include "simplemd/configurations/ExternalForceConfiguration.h"
#include "simplemd/configurations/MPIConfiguration.h"
#include "simplemd/configurations/MoleculeConfiguration.h"
#include "simplemd/configurations/ProfilePlotterConfiguration.h"
#include "simplemd/configurations/RDFConfiguration.h"
#include "simplemd/configurations/SimulationConfiguration.h"
#include "simplemd/configurations/VTKConfiguration.h"
#include "tarch/configuration/Configuration.h"
#include <iostream>
#include <vector>

namespace simplemd {
namespace configurations {
class MolecularDynamicsConfiguration;
}
} // namespace simplemd

/** configuration for MD simulations. Bundles all other configuration tags in
 * one main tag.
 *  @author Philipp Neumann
 */
class simplemd::configurations::MolecularDynamicsConfiguration : public tarch::configuration::Configuration {

public:
  MolecularDynamicsConfiguration() {}
  virtual ~MolecularDynamicsConfiguration() {}

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

  const simplemd::configurations::DomainConfiguration& getDomainConfiguration() const { return _domainConfiguration; }

  simplemd::configurations::DomainConfiguration& getDomainConfigurationNonConst() { return _domainConfiguration; }

  const simplemd::configurations::MoleculeConfiguration& getMoleculeConfiguration() const { return _moleculeConfiguration; }

  const simplemd::configurations::MPIConfiguration& getMPIConfiguration() const { return _mpiConfiguration; }

  const simplemd::configurations::VTKConfiguration& getVTKConfiguration() const { return _vtkConfiguration; }

  const simplemd::configurations::SimulationConfiguration& getSimulationConfiguration() const { return _simulationConfiguration; }

  const simplemd::configurations::RDFConfiguration& getRDFConfiguration() const { return _rdfConfiguration; }

  const simplemd::configurations::CheckpointConfiguration& getCheckpointConfiguration() const { return _checkpointConfiguration; }

  const std::vector<simplemd::configurations::ProfilePlotterConfiguration>& getProfilePlotterConfigurations() const { return _profilePlotterConfigurations; }

  const std::vector<simplemd::configurations::ExternalForceConfiguration>& getExternalForceConfigurations() const { return _externalForceConfigurations; }

private:
  simplemd::configurations::DomainConfiguration _domainConfiguration;
  simplemd::configurations::MoleculeConfiguration _moleculeConfiguration;
  simplemd::configurations::MPIConfiguration _mpiConfiguration;
  simplemd::configurations::VTKConfiguration _vtkConfiguration;
  simplemd::configurations::SimulationConfiguration _simulationConfiguration;
  simplemd::configurations::RDFConfiguration _rdfConfiguration;
  simplemd::configurations::CheckpointConfiguration _checkpointConfiguration;
  std::vector<simplemd::configurations::ProfilePlotterConfiguration> _profilePlotterConfigurations;
  std::vector<simplemd::configurations::ExternalForceConfiguration> _externalForceConfigurations;
};

#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_MOLECULARDYNAMICSCONFIGURATION_H_
