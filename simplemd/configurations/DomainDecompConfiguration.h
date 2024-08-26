// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#ifndef _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINDECOMPCONFIGURATION_H_
#define _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINDECOMPCONFIGURATION_H_

#include "simplemd/MolecularDynamicsDefinitions.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/la/Vector.h"
#include <array>
#include <iostream>
#include <vector>

namespace simplemd {
namespace configurations {
class DomainDecompConfiguration;
}
} // namespace simplemd

/** input configuration for specifying domain decomposition type
 *  @author Amartya Das Sharma
 */

class simplemd::configurations::DomainDecompConfiguration : public tarch::configuration::Configuration {
public:
  DomainDecompConfiguration();
  virtual ~DomainDecompConfiguration() {}

  /** enum containing the type of domain decompositions currently supported
   * DEFAULT : default behaviour, regularly spaced grid
   * STATIC_IRREG_RECT_GRID : rectilinear grid with weights given to specify subdomain sizes
   * E.g. x values 1,2,1 define an x axis with subdomain lengths in the 1:2:1 ratio
   */
  enum class DecompositionType { DEFAULT, STATIC_IRREG_RECT_GRID };

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
   * If a tag isn't optional and parseSubtag() was not called (first case)
   */
  bool isValid() const;

  /** getters for all parsed and computed quantities */
  const tarch::la::Vector<MD_DIM, std::vector<unsigned int>>& getSubdomainWeights() const { return _subdomainWeights; }
  DecompositionType getDecompType() const { return _decompType; }

private:
  static const std::string DECOMP_TYPE;
  static const std::string DEFAULT_DECOMP;
  static const std::string STATIC_IRREG_RECT_GRID;
  static const std::string AXES[MD_DIM];

  /** helper method to parse input string */
  std::vector<unsigned int> getWeightsFromString(std::string weights);

  /** type of decomposition */
  DecompositionType _decompType;

  /** weights parsed from the XML */
  tarch::la::Vector<MD_DIM, std::vector<unsigned int>> _subdomainWeights;

  /** whether the decomposition-type tag was present or not (since the tag is optional) */
  bool _isDefined;

  /** validity of the config */
  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_CONFIGURATIONS_DOMAINDECOMPCONFIGURATION_H_