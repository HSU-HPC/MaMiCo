// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_VTKCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_VTKCONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace simplemd {
namespace configurations {
class VTKConfiguration;
}
} // namespace simplemd

/** configuration input for VTK output.
 *  @author Philipp Neumann
 */
class simplemd::configurations::VTKConfiguration : public tarch::configuration::Configuration {
public:
  VTKConfiguration();
  virtual ~VTKConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement* node);

  /**
   * Return name of xml tag that is associated to the configuration.
   */
  virtual std::string getTag() const;

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
  const std::string& getFilename() const { return _filename; }
  const unsigned int& getWriteEveryTimestep() const { return _writeEveryTimestep; }

private:
  static const std::string FILENAME;
  static const std::string WRITE_EVERY_TIMESTEP;

  /** name stem of the corresponding vtk file */
  std::string _filename;

  /** number of timesteps per vtk output. If this value is zero, no output is
   * written at all */
  unsigned int _writeEveryTimestep;

  bool _isValid;
};
#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_VTKCONFIGURATION_H_
