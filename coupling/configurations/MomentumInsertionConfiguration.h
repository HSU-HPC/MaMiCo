// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MOMENTUMINSERTIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MOMENTUMINSERTIONCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>
#include "coupling/AdditiveMomentumInsertion.h"
#include "coupling/SetGivenVelocity4MomentumInsertion.h"
#include "coupling/VelocityGradientRelaxation.h"
#include "coupling/NoMomentumInsertion.h"
#include "coupling/NieVelocityImposition.h"

namespace coupling {
  namespace configurations {
    class MomentumInsertionConfiguration;
  }
}

// forward declaration for testing
class NieTest;

/** reads momentum insertion tag.
 *  @author Philipp Neumann
 */
class coupling::configurations::MomentumInsertionConfiguration:
public tarch::configuration::Configuration {
  public:
    // for testing: give access to variables; will only be used in the read-sense
    friend class ::NieTest;

    enum MomentumInsertionType{
      ADDITIVE_MOMENTUM_INSERTION=0,DIRECT_VELOCITY_INSERTION=1,VELOCITY_GRADIENT_RELAXATION=2,NO_INSERTION=3,
      VELOCITY_GRADIENT_RELAXATION_TOPONLY=4,NIE_VELOCITY_IMPOSITION=5
    };

    MomentumInsertionConfiguration(): _insertionType(ADDITIVE_MOMENTUM_INSERTION),_isValid(true){}

    virtual ~MomentumInsertionConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node ){
      std::string value;
      tarch::configuration::ParseConfiguration::readStringMandatory(value,node,"type");
      if (value=="additive-momentum-insertion"){
        _insertionType = ADDITIVE_MOMENTUM_INSERTION;
      } else if (value=="direct-velocity-insertion"){
        _insertionType = DIRECT_VELOCITY_INSERTION;
      } else if (value=="velocity-gradient-relaxation"){
        _insertionType = VELOCITY_GRADIENT_RELAXATION;
      } else if (value=="none"){
        _insertionType = NO_INSERTION;
      } else if (value=="velocity-gradient-relaxation-top-only"){
        _insertionType = VELOCITY_GRADIENT_RELAXATION_TOPONLY;
      } else if (value=="nie-velocity-imposition"){
        _insertionType = NIE_VELOCITY_IMPOSITION;
      } else {
        std::cout << "ERROR coupling::MomentumInsertionConfiguration: Wrong insertion type!" << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }

      if (_insertionType==VELOCITY_GRADIENT_RELAXATION){
        tarch::configuration::ParseConfiguration::readDoubleMandatory(_velocityRelaxationFactor,node,"velocity-relaxation-factor");
        if ((_velocityRelaxationFactor<=0.0) || (_velocityRelaxationFactor > 1.0) ){
          std::cout << "ERROR coupling::MomentumInsertionConfiguration: velocity-relaxation-factor=" << _velocityRelaxationFactor << "!";
          std::cout << " It must be in the range (0.0,1.0]!" << std::endl;
          _isValid = false;
          exit(EXIT_FAILURE);
        }
      } else if (_insertionType==VELOCITY_GRADIENT_RELAXATION_TOPONLY){
        tarch::configuration::ParseConfiguration::readDoubleMandatory(_velocityRelaxationFactor,node,"velocity-relaxation-factor");
        if ((_velocityRelaxationFactor<=0.0) || (_velocityRelaxationFactor > 1.0) ){
          std::cout << "ERROR coupling::MomentumInsertionConfiguration: velocity-relaxation-factor=" << _velocityRelaxationFactor << "!";
          std::cout << " It must be in the range (0.0,1.0]!" << std::endl;
          _isValid = false;
          exit(EXIT_FAILURE);
        }
      } else if (_insertionType==NIE_VELOCITY_IMPOSITION){
        int buf;
        tarch::configuration::ParseConfiguration::readIntMandatory(buf,node,"innermost-overlap-layer");
        if (buf <=0){std::cout << "ERROR coupling::MomentumInsertionConfiguration: innermost-overlap-layer=" << buf <<"!" << std::endl; _isValid=false; exit(EXIT_FAILURE);}
        _innerOverlap = (unsigned int) buf;
        tarch::configuration::ParseConfiguration::readIntMandatory(buf,node,"outermost-overlap-layer");     
        if ( ((unsigned int)buf > _innerOverlap) || (buf <= 0) ){std::cout << "ERROR coupling::MomentumInsertionConfiguration: outermost-overlap-layer=" << buf <<"!" << std::endl; _isValid=false; exit(EXIT_FAILURE);}
        _outerOverlap = (unsigned int) buf;

      }
    }

    /**
     * Return name of xml tag that is associated to the configuration.
     */
    std::string getTag() const {return "momentum-insertion";}

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

    /** getters */
    const MomentumInsertionType& getMomentumInsertionType() const {return _insertionType;}

    template<class LinkedCell,unsigned int dim>
    MomentumInsertion<LinkedCell,dim>* interpreteConfiguration(
      coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface,
      const coupling::IndexConversion<dim>& indexConversion,
      const coupling::datastructures::MacroscopicCellWithLinkedCells<LinkedCell,dim> * const macroscopicCells,
      unsigned int numberMDTimestepsPerCouplingCycle
    ) const {
      if (_insertionType == ADDITIVE_MOMENTUM_INSERTION){
        return new coupling::AdditiveMomentumInsertion<LinkedCell,dim>(mdSolverInterface,numberMDTimestepsPerCouplingCycle);
      } else if (_insertionType == DIRECT_VELOCITY_INSERTION){
        return new coupling::SetGivenVelocity4MomentumInsertion<LinkedCell,dim>(mdSolverInterface);
      } else if (_insertionType == VELOCITY_GRADIENT_RELAXATION){
        return new coupling::VelocityGradientRelaxation<LinkedCell,dim>(_velocityRelaxationFactor,mdSolverInterface,indexConversion,macroscopicCells);
      } else if (_insertionType == NO_INSERTION){
        return new coupling::NoMomentumInsertion<LinkedCell,dim>(mdSolverInterface);
      } else if (_insertionType == VELOCITY_GRADIENT_RELAXATION_TOPONLY){
        return new coupling::VelocityGradientRelaxationTopOnly<LinkedCell,dim>(_velocityRelaxationFactor,mdSolverInterface,indexConversion,macroscopicCells);
      } else if (_insertionType == NIE_VELOCITY_IMPOSITION){
        return new coupling::NieVelocityImposition<LinkedCell,dim>(mdSolverInterface,indexConversion,_outerOverlap,_innerOverlap);
      }
      return NULL;
    }

  private:
    MomentumInsertionType _insertionType;
    double _velocityRelaxationFactor;    // required by velocity relaxation schemes
    unsigned int _innerOverlap;          // innermost layer of overlap cells (required by nie velocity imposition)
    unsigned int _outerOverlap;          // outermost layer of overlap cells (required by nie velocity imposition)
    bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_MOMENTUMINSERTIONCONFIGURATION_H_
