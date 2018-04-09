#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <iostream>
#include "coupling/noisereduction/NoiseReduction.h"
#include "coupling/noisereduction/IdentityTransform.h"

namespace coupling {
  namespace configurations {
    class NoiseReductionConfiguration;
  }
}

/** noise reduction configuration, i.e. smoothing algorithm to filter MD fluctuations
 *  @author Piet Jarmatz
 */
class coupling::configurations::NoiseReductionConfiguration:
public tarch::configuration::Configuration {
  public:
    enum NoiseReductionType{
      IdentityTransform=0, GaussianFilter=1, POD=2, AnisotropicDiffusion=3
    };

    NoiseReductionConfiguration(): _type(IdentityTransform), _isValid(true){}

    virtual ~NoiseReductionConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node ){
      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value,node,"type");
      if (value=="none"){
        _type = IdentityTransform;
      } else if (value=="gaussian-filter"){
        _type = GaussianFilter;
      }  else if (value=="POD"){
        _type = POD;
      } else if (value=="anisotropic-diffusion"){
        _type = AnisotropicDiffusion;
      } else {
        std::cout << "ERROR coupling::NoiseReductionConfiguration: Wrong noise filter type!" << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
    }

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "noise-reduction";}

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


    template<class LinkedCell, unsigned int dim>
    coupling::noisereduction::NoiseReduction<LinkedCell,dim>* interpreteConfiguration() const {
      if (_type == IdentityTransform){
        return new coupling::noisereduction::IdentityTransform<LinkedCell,dim>();
      } else if (_type == GaussianFilter){
        std::cout << "ERROR coupling::NoiseReductionConfiguration: not implemented" << std::endl;
        exit(EXIT_FAILURE);
      } else if (_type == POD){
        std::cout << "ERROR coupling::NoiseReductionConfiguration: not implemented" << std::endl;
        exit(EXIT_FAILURE);
      } else if (_type == AnisotropicDiffusion){
        std::cout << "ERROR coupling::NoiseReductionConfiguration: not implemented" << std::endl;
        exit(EXIT_FAILURE);
      } else {
        return NULL;
      }
    }

    NoiseReductionType getNoiseReductionType() const {return _type;}

  private:
    NoiseReductionType _type;

    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
