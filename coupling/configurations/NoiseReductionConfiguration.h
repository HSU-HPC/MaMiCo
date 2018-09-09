#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include <iostream>
#include "coupling/noisereduction/NoiseReduction.h"
#include "coupling/noisereduction/IdentityTransform.h"
#include "coupling/noisereduction/POD.h"

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

        int buf = -1;

        tarch::configuration::ParseConfiguration::readIntMandatory(buf,node,"time-window-size");
        if (buf <= 2){
          std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
          std::cout << "parseSubtag(): " << "time-window-size" << " smaller than or equal two!" << std::endl;
          exit(EXIT_FAILURE);
        } else {
          _tws = buf;
        }

        tarch::configuration::ParseConfiguration::readIntMandatory(buf,node,"kmax");
        if (buf <= 0){
          std::cout << "ERROR coupling::configurations::ParticleInsertionConfiguration::";
          std::cout << "parseSubtag(): " << "kmax" << " smaller or equal zero!" << std::endl;
          exit(EXIT_FAILURE);
        } else {
          _kmax = buf;
        }

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

    // tws_param can be used to override XML configuration
    template<unsigned int dim>
    coupling::noisereduction::NoiseReduction<dim>* interpreteConfiguration(
      const coupling::IndexConversion<dim> &indexConversion,
      const tarch::utils::MultiMDService<dim>& multiMDService,
      int tws_param = 0
    ) const {
      if (_type == IdentityTransform){
        return new coupling::noisereduction::IdentityTransform<dim>(indexConversion, multiMDService);
      } else if (_type == GaussianFilter){
        std::cout << "ERROR coupling::NoiseReductionConfiguration: not implemented" << std::endl;
        exit(EXIT_FAILURE);
      } else if (_type == POD){
        return new coupling::noisereduction::POD<dim>(indexConversion, multiMDService, tws_param==0?_tws:tws_param, _kmax);
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
    int _tws;
    int _kmax;

    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_NOISEREDUCTIONCONFIGURATION_H_
