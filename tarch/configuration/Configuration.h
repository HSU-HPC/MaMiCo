// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_CONFIGURATION_CONFIGURATION_H_
#define _TARCH_CONFIGURATION_CONFIGURATION_H_

#include "tarch/tinyxml2/tinyxml2.h"
#include <string>

namespace tarch {
/**	namespace tarch */
namespace configuration {
/**	namespace configuration */
class Configuration;
}
}

/**	interface for configuration using tinyxml2
 *  @author Philipp Neumann
 */
class tarch::configuration::Configuration {
public:
  /**
   * Destructor
   */
  virtual ~Configuration() {}

  /** Parses a XML-element.
	 * @param node a pointer to the xml-element from tinyxml2 class
     */
  virtual void parseSubtag(tinyxml2::XMLElement *node) = 0;

  /**
   * @return the name of the current xml-element
   */
  virtual std::string getTag() const = 0;

  /**
   *  @returns true if the configuration is valid; false otherwise.
   *  This method is basically for debugging purposes and detection of invalid
   * input.
   */
  virtual bool isValid() const = 0;
};

#endif
