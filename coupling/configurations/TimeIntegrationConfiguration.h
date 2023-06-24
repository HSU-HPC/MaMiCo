// Copyright (C) 2023 Helmut Schmidt University
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#pragma once

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <iostream>

namespace coupling {
namespace configurations {
class TimeIntegrationConfiguration;
}
} // namespace coupling

/** 
 *  Reads time integration configuration. Use number-subdomains="1" to disable PinT. 
 *  Derived from tarch::configuration::Configuration
 *	@brief reads time integration configuration
 *  @author Piet Jarmatz
 */
class coupling::configurations::TimeIntegrationConfiguration : public tarch::configuration::Configuration {
public:
  /** Constructor, initializes the class  */
  TimeIntegrationConfiguration() : _pint_domains(1), _pint_iterations(1), _auto_iteration(false), _isValid(true) {}

  /** Destructor */
  virtual ~TimeIntegrationConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement* node) {
    int buf;
    tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, "number-subdomains");
    if (buf <= 0) {
      std::cout << "ERROR coupling::TimeIntegrationConfiguration: "
                   "number-subdomains = "
                << buf << "!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
    _pint_domains = buf;

    if(isPinTEnabled()){
      #if (COUPLING_MD_PARALLEL != COUPLING_MD_YES)
        std::cout << "ERROR coupling::TimeIntegrationConfiguration: PinT is enabled but COUPLING_MD_PARALLEL disabled" << std::endl;
        std::cout << "Disable PinT in XML config, or enable BUILD_WITH_MPI in cmake." << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      #endif

      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value, node, "number-iterations");
      if (value == "auto"){
        _auto_iteration = true;
        _pint_iterations = 0;
      }
      else{
        _auto_iteration = false;
        tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, "number-iterations");
        if (buf <= 0) {
          std::cout << "ERROR TimeIntegrationConfiguration::number-iterations too small!" << std::endl;
          _isValid = false;
          exit(EXIT_FAILURE);
        }
        if (buf > _pint_domains) {
          std::cout << "ERROR TimeIntegrationConfiguration::number-iterations too large!" << std::endl;
          _isValid = false;
          exit(EXIT_FAILURE);
        }
        _pint_iterations = buf;
      }  
    }
  }

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "time-integration"; }

  /** checks if the configuration is valid. 
   * 	@return _isValid
   */
  bool isValid() const { return _isValid; }

  int getPintDomains() const { return _pint_domains; }
  bool isPinTEnabled() const { return _pint_domains > 1; }
  int getPintIterations() const { return _pint_iterations; }
  bool isAutoIteration() const { return _auto_iteration; }

private:
  int _pint_domains;
  int _pint_iterations;
  bool _auto_iteration;

  bool _isValid;
};
