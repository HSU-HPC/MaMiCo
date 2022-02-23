// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_MPICONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_MPICONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace simplemd {
namespace configurations {
class MPIConfiguration;
}
} // namespace simplemd

/** configuration for distributed MPI MD simulations. Currently only contains
 * the splitting of the processes for the domain.
 *  @author Philipp Neumann
 */
class simplemd::configurations::MPIConfiguration
    : public tarch::configuration::Configuration {
public:
  MPIConfiguration();
  virtual ~MPIConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement *node);

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

  /** getters */
  const tarch::la::Vector<MD_DIM, unsigned int> &getNumberOfProcesses() const {
    return _numberOfProcesses;
  }

private:
  static const std::string NUMBER_OF_PROCESSES;

  tarch::la::Vector<MD_DIM, unsigned int> _numberOfProcesses;
  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_MPICONFIGURATION_H_
