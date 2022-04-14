// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_

#include "coupling/paralleltopology/ParallelTopologyFactory.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace coupling {
namespace configurations {
class ParallelTopologyConfiguration;
}
} // namespace coupling

/** reads parallel topology configuration. XYZ and ZYX are supported.. Derive
 *from the class tarch::configuration::Configuration
 *	@brief reads parallel topology configuration. XYZ and ZYX are supported.
 *  @author Philipp Neumann
 */
class coupling::configurations::ParallelTopologyConfiguration : public tarch::configuration::Configuration {
public:
  /** Constructor, initializes the class  */
  ParallelTopologyConfiguration() : _type(coupling::paralleltopology::UNDEFINED), _isValid(true) {}

  /** Destructor */
  virtual ~ParallelTopologyConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement* node) {
    std::string value;
    tarch::configuration::ParseConfiguration::readStringMandatory(value, node, "type");
    if (value == "xyz") {
      _type = coupling::paralleltopology::XYZ;
    } else if (value == "zyx") {
      _type = coupling::paralleltopology::ZYX;
    } else {
      std::cout << "ERROR coupling::ParallelTopologyConfiguration: Wrong type!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
  }

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "parallel-topology"; }

  /** checks if the configuration is valid. This operation usually fails, if
   *e.g.
   *	1. parseSubtag() hasn't been called, i.e. configuration has not been
   *used, or
   *  2. parseSubtag() failed due to a wrong file.
   * 	@return _isValid
   */
  bool isValid() const { return _isValid; }

  /** Returns the parallel topology type.
   * 	@return _type
   */
  coupling::paralleltopology::ParallelTopologyType getParallelTopologyType() const { return _type; }

private:
  coupling::paralleltopology::ParallelTopologyType _type;

  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_
