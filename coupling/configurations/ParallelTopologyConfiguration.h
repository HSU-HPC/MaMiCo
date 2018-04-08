// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>
#include "coupling/paralleltopology/ParallelTopologyFactory.h"

namespace coupling {
  namespace configurations {
    class ParallelTopologyConfiguration;
  }
}


/** reads parallel topology configuration. XYZ and ZYX are supported.
 *  @author Philipp Neumann
 */
class coupling::configurations::ParallelTopologyConfiguration:
public tarch::configuration::Configuration {
  public:

    ParallelTopologyConfiguration(): _type(coupling::paralleltopology::UNDEFINED),_isValid(true){}

    virtual ~ParallelTopologyConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node ){
      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value,node,"type");
      if (value=="xyz"){
        _type = coupling::paralleltopology::XYZ;
      } else if (value=="zyx"){
        _type = coupling::paralleltopology::ZYX;
      } else {
        std::cout << "ERROR coupling::ParallelTopologyConfiguration: Wrong type!" << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
    }

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "parallel-topology";}

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

    coupling::paralleltopology::ParallelTopologyType getParallelTopologyType() const { return _type;}

  private:
    coupling::paralleltopology::ParallelTopologyType _type;

    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARALLELTOPOLOGYCONFIGURATION_H_
