// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARTICLEINSERTIONCONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARTICLEINSERTIONCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include "coupling/interface/MDSolverInterface.h"
#include "coupling/UsherParticleInsertion.h"
#include "coupling/NoParticleInsertion.h"
#include <iostream>

namespace coupling {
  namespace configurations {
    class ParticleInsertionConfiguration;
  }
}


/** configuration for particle insertion algorithm (e.g.: USHER).
 *  @author Philipp Neumann
 */
class coupling::configurations::ParticleInsertionConfiguration: public tarch::configuration::Configuration {
  public:
    ParticleInsertionConfiguration():
    _insertDeleteMassEveryTimestep(1),
    _rSigmaCoeff(0.0),
    _meanPotentialEnergyFactor(0.0),
    _uOverlapCoeff(0.0),
    _stepRefCoeff(0.0),
    _iterMax(0),
    _restartMax(0),
    _tolerance(0.0),
    _offsetFromOuterBoundary(0.0){}

    virtual ~ParticleInsertionConfiguration(){}

    void parseSubtag( tinyxml2::XMLElement* node );

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
     * If a tag ain't optional and parseSubtag() was not called (first case)
     */
    bool isValid() const;

    template<class LinkedCell, unsigned int dim>
    coupling::ParticleInsertion<LinkedCell,dim>* interpreteConfiguration(coupling::interface::MDSolverInterface<LinkedCell,dim> * const mdSolverInterface) const {
      if (_particleInsertionType==USHER){
        return new coupling::UsherParticleInsertion<LinkedCell,dim>(
          _insertDeleteMassEveryTimestep,_rSigmaCoeff, _meanPotentialEnergyFactor, _uOverlapCoeff, _stepRefCoeff, _iterMax, _restartMax, _tolerance, _offsetFromOuterBoundary,
          mdSolverInterface);
      } else if (_particleInsertionType==NO_INSERTION){
        return new coupling::NoParticleInsertion<LinkedCell,dim>();
      }

      // this case should never appear
      return NULL;
    }

    enum ParticleInsertionType{USHER=0,NO_INSERTION=1};
    ParticleInsertionType getParticleInsertionType() const { return _particleInsertionType; }

  private:
    static const std::string INSERT_DELETE_MASS_EVERY_TIMESTEP;
    static const std::string RSIGMA_COEFF;
    static const std::string MEAN_POTENTIAL_ENERGY_FACTOR;
    static const std::string UOVERLAP_COEFF;
    static const std::string STEPREF_COEFF;
    static const std::string ITER_MAX;
    static const std::string RESTART_MAX;
    static const std::string TOLERANCE;
    static const std::string OFFSET_FROM_OUTER_BOUNDARY;

    /** used to trigger mass insertion/deletion actions only every X time steps. For NoParticleInsertion, a value of 1 is used since this is irrelevant in this case.
     *  If the value is not specified in the configuration, a default value of 1 is used. */
    unsigned int _insertDeleteMassEveryTimestep;

    /** r_sigma value. In the USHER paper, this value resembles the respective parameter for choosing the stepsize
     *  when two particles overlap. In the paper, it is stated that it should be close to unity (in sigma units; the
     *  parameter here is also used in sigma units of the respective molecules) and is chosen to be 0.9. If the parameter
     *  is not defined in the config, _rSigmaCoeff is also set to 0.9.
     */
    double _rSigmaCoeff;
    double _meanPotentialEnergyFactor;

    /** overlap energy. If not defined in the config, this value is set to 10^4. It is typically given
     *  in units of epsilon (energy parameter of the LJ fluid).
     */
    double _uOverlapCoeff;

    /** coefficient determining a reference step size when stepping towards the right energy level.
     *  From the USHER paper (Usher - An algorithm for particle insertion in dense fluids), this
     *  coefficient is optimally chosen to be 0.1 (in sigma units; this configuration also determines
     *  this parameter in sigma units). Thus, the parameter is set to 0.1, if it is not defined
     *  within the configuration.
     */
    double _stepRefCoeff;

    /** maximum number of USHER iteration steps within one insertion cycle */
    unsigned int _iterMax;

    /** maximum number of USHER restarts, that is USHER insertion cycles */
    unsigned int _restartMax;

    /** tolerance for USHER stopping criterion. If this is not specified, we
     *  set it to _meanPotentialEnergyFactor*2.0.
     */
    double _tolerance;

    /** enables valid particle positions only if the particle is at least at this distance from an outer boundary.
     *  This is useful since open boundary forces are not considered in the current USHER implementation. Inserting a particle
     *  close to an outer open boundary may thus result in very strong forces and, thus, yield instabilities.
     */
    double _offsetFromOuterBoundary;

    /** type of particle insertion algorithm */
    ParticleInsertionType _particleInsertionType;
};
#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_PARTICLEINSERTIONCONFIGURATION_H_
