// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_THERMOSTATCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_THERMOSTATCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

namespace coupling {
  namespace configurations {
    class ThermostatConfiguration;
  }
}

/** reads the configuartion for the domain of the thermostat.
 *  @author Helene Wittenberg
 */
class coupling::configurations::ThermostatConfiguration:
public tarch::configuration::Configuration {
  public:
    enum ThermostatRegion{outer, overlap, all, none};
    ThermostatConfiguration(): _type{none},_isValid(true){}

    virtual ~ThermostatConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node ){
      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value,node,"type");
      if (value=="outer"){
        _type = outer;
      } else if (value=="overlap"){
        _type = overlap;
      } else if (value=="all"){
        _type = all;
      } else {
        std::cout << "ERROR coupling::ThermostatConfiguration: Wrong type!" << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
    }

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "thermostat";}

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
    bool isValid() const { return _isValid;}

    ThermostatRegion getThermostatRegionType() const { return _type;}

  private:
    ThermostatRegion _type;

    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_THERMOSTATCONFIGURATION_H_
