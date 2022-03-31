#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_

#include "coupling/noisereduction/IdentityTransform.h"
#include "coupling/noisereduction/NLM.h"
#include "coupling/noisereduction/NoiseReduction.h"
#include "coupling/noisereduction/POD.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <iostream>

namespace coupling {
namespace configurations {
class NoiseReductionConfiguration;
}
} // namespace coupling

/** noise reduction configuration, i.e. smoothing algorithm to filter MD
 * fluctuations. Derive from the class tarch::configuration::Configuration
 *	@brief noise reduction configuration, i.e. smoothing algorithm to filter
 *MD fluctuations
 *  @author Piet Jarmatz
 *	@ Piet could you please take a look on this class
 */
class coupling::configurations::NoiseReductionConfiguration : public tarch::configuration::Configuration {
public:
  /** noise reduction types that are implemented.
   *	@enum NoiseReductionType
   */
  enum NoiseReductionType {
    IdentityTransform = 0 /**< no filtering*/,
    GaussianFilter = 1 /**< Gaussian filter*/,
    POD = 2 /**< Proper orthogonal decomposition filter*/,
    NLM = 3 /**< non local means filter*/
  };

  /** Constructor, initializes the class  */
  NoiseReductionConfiguration() : _type(IdentityTransform), _isValid(true) {}

  /** Destructor */
  virtual ~NoiseReductionConfiguration() {}

  void parseSubtag(tinyxml2::XMLElement *node) {
    std::string value;
    tarch::configuration::ParseConfiguration::readStringMandatory(value, node, "type");
    if (value == "none") {
      _type = IdentityTransform;
    } else if (value == "gaussian-filter") {
      _type = GaussianFilter;
    } else if (value == "POD") {
      _type = POD;

      int buf = -1;

      tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, "time-window-size");
      if (buf <= 2) {
        std::cout << "ERROR "
                     "coupling::configurations::ParticleInsertionConfiguration::";
        std::cout << "parseSubtag(): "
                  << "time-window-size"
                  << " smaller than or equal two!" << std::endl;
        exit(EXIT_FAILURE);
      } else {
        _tws = buf;
      }

      tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, "kmax");
      if (buf <= 0) {
        std::cout << "ERROR "
                     "coupling::configurations::ParticleInsertionConfiguration::";
        std::cout << "parseSubtag(): "
                  << "kmax"
                  << " smaller or equal zero!" << std::endl;
        exit(EXIT_FAILURE);
      } else {
        _kmax = buf;
      }

    } else if (value == "NLM") {
      _type = NLM;

      int buf = -1;

      tarch::configuration::ParseConfiguration::readIntMandatory(buf, node, "time-window-size");
      if (buf <= 2) {
        std::cout << "ERROR "
                     "coupling::configurations::ParticleInsertionConfiguration::";
        std::cout << "parseSubtag(): "
                  << "time-window-size"
                  << " smaller than or equal two!" << std::endl;
        exit(EXIT_FAILURE);
      } else {
        _tws = buf;
      }

    } else {
      std::cout << "ERROR coupling::NoiseReductionConfiguration: Wrong noise "
                   "filter type!"
                << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
  }

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "noise-reduction"; }

  /** checks if the configuration is valid. This operation usually fails, if
e.g.
         *	1. parseSubtag() hasn't been called, i.e. configuration has not
been used, or
     *  2. parseSubtag() failed due to a wrong file.
         * 	@return _isValid
     */
  bool isValid() const { return _isValid; }

  /** Returns noise reduction configuration.
   * 	@tparam dim Number of dimensions; it can be 1, 2 or 3
   * 	@param indexConversion
   * 	@param multiMDService
   * 	@param tws_param =0 by default
   * 	@return noise reduction config
   *	@note tws_param can be used to override XML configuration
   */
  template <unsigned int dim>
  coupling::noisereduction::NoiseReduction<dim> *interpreteConfiguration(const coupling::IndexConversion<dim> &indexConversion,
                                                                         const tarch::utils::MultiMDService<dim> &multiMDService, int tws_param = 0) const {
    if (_type == IdentityTransform) {
      return new coupling::noisereduction::IdentityTransform<dim>(indexConversion, multiMDService);
    } else if (_type == GaussianFilter) {
      std::cout << "ERROR coupling::NoiseReductionConfiguration: not implemented" << std::endl;
      exit(EXIT_FAILURE);
    } else if (_type == POD) {
      return new coupling::noisereduction::POD<dim>(indexConversion, multiMDService, tws_param == 0 ? _tws : tws_param, _kmax);
    } else if (_type == NLM) {
      return new coupling::noisereduction::NLM<dim>(indexConversion, multiMDService, _tws);
    } else {
      return NULL;
    }
  }

  /** Returns the noise reduction type.
   * 	@return _type
   */
  NoiseReductionType getNoiseReductionType() const { return _type; }

private:
  NoiseReductionType _type;
  int _tws;
  int _kmax;

  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
