// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_STRATEGYCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_STRATEGYCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>
#include "coupling/transferstrategies/TransferStrategy.h"
#include "coupling/transferstrategies/DifferenceTransferStrategy.h"
#include "coupling/transferstrategies/DirectTransferStrategy.h"
#include "coupling/transferstrategies/TransferStrategy4SchwarzCoupling.h"
#include "coupling/transferstrategies/TransferStrategy4NieCoupling.h"
#include "coupling/transferstrategies/AveragingTransferStrategy.h"

namespace coupling {
  namespace configurations {
    template<unsigned int dim>
    class TransferStrategyConfiguration;
  }
}


/** transfer strategy configuration, i.e. algorithm/combin. of quantity transfer steps and quantity interpretation (e.g. momentum vs. velocity).
 *  @author Philipp Neumann
 */
template<unsigned int dim>
class coupling::configurations::TransferStrategyConfiguration:
public tarch::configuration::Configuration {
  public:
    enum StrategyType{
      DirectTransferStrategy=0,DifferenceTransferStrategy=1,TransferStrategy4FluxCoupling=2,TransferStrategy4SchwarzCoupling=3,AveragingTransferStrategy=4,
      TransferStrategy4NieCoupling=5
    };

    TransferStrategyConfiguration(): _type(DirectTransferStrategy),_massFluxBoundary(false),_isValid(true){}

    virtual ~TransferStrategyConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node ){
      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value,node,"type");
      if (value=="direct-transfer"){
        _type = DirectTransferStrategy;
      } else if (value=="difference-transfer"){
        _type = DifferenceTransferStrategy;
      }  else if (value=="schwarz-transfer"){
        _type = TransferStrategy4SchwarzCoupling;
      } else if (value=="nie-transfer"){
        _type = TransferStrategy4NieCoupling;
      } else if (value=="averaging"){
        _type = AveragingTransferStrategy;
      } else {
        std::cout << "ERROR coupling::TransferStrategyConfiguration: Wrong insertion type!" << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }


      if (_type==TransferStrategy4NieCoupling){
        const std::string boundaries[6] = {"mass-flux-west","mass-flux-east","mass-flux-south","mass-flux-north","mass-flux-bottom","mass-flux-top"};
        for (unsigned int d = 0; d < 2*dim; d++){
          tarch::configuration::ParseConfiguration::readBoolMandatory(_massFluxBoundary[d],node,boundaries[d]);
        }
        // by default: no shifting of time interval
        _shiftTimestep=0.0;
        tarch::configuration::ParseConfiguration::readDoubleOptional(_shiftTimestep,node,"shift-by-timesteps");
        if (_shiftTimestep<0.0 || _shiftTimestep>1.0){
          std::cout << "Warning coupling::configurations::TransferStrategyConfiguration: shift-by-timesteps=" << _shiftTimestep << "; typical values range between 0 and 1!" << std::endl;
        }
      }
    }

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "transfer-strategy";}

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


    template<class LinkedCell>
    coupling::transferstrategies::TransferStrategy<LinkedCell,dim>* interpreteConfiguration(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim> &indexConversion,
      unsigned int numberOfMDTimesteps
    ) const {
      if (_type == DirectTransferStrategy){
        return new coupling::transferstrategies::DirectTransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion);
      } else if (_type == DifferenceTransferStrategy){
        return new coupling::transferstrategies::DifferenceTransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion,numberOfMDTimesteps);
      } else if (_type == TransferStrategy4SchwarzCoupling){
        return new coupling::transferstrategies::TransferStrategy4SchwarzCoupling<LinkedCell,dim>(mdSolverInterface,indexConversion,numberOfMDTimesteps);
      } else if (_type == TransferStrategy4NieCoupling){
        return new coupling::transferstrategies::TransferStrategy4NieCoupling<LinkedCell,dim>(mdSolverInterface,indexConversion,numberOfMDTimesteps,_shiftTimestep,_massFluxBoundary);
      } else if (_type == AveragingTransferStrategy){
        return new coupling::transferstrategies::AveragingTransferStrategy<LinkedCell,dim>(mdSolverInterface,indexConversion);
      } else {
        return NULL;
      }
    }

    StrategyType getStrategyType() const {return _type;}

  private:
    StrategyType _type;
    tarch::la::Vector<2*dim,bool> _massFluxBoundary; // true in each component, if one of the 2*dim boundaries allows for mass flux
    double _shiftTimestep; // used for Nie coupling: time interval by which the evaluation of the continuum flow field should be shifted. See also TransferStrategy4NieCoupling for more details.

    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_STRATEGYCONFIGURATION_H_
