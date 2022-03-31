// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_BOUNDARYFORCECONFIGURATION_H_
#define _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_BOUNDARYFORCECONFIGURATION_H_

#include "coupling/NoBoundaryForce.h"
#include "coupling/ZhouBoundaryForceController.h"
#include "coupling/interface/MDSolverInterface.h"
#include "tarch/configuration/Configuration.h"
#include "tarch/configuration/ParseConfiguration.h"
#include "tarch/la/Vector.h"
#include <iostream>

namespace coupling {
namespace configurations {
template <unsigned int dim> class BoundaryForceConfiguration;
}
} // namespace coupling

/** boundary force configuration
 *	@brief reads boundary force tag
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim> class coupling::configurations::BoundaryForceConfiguration : public tarch::configuration::Configuration {
public:
  /** boundary force types that are implemented.
   *	@enum BoundaryForceType
   */
  enum BoundaryForceType {
    ZHOU = 0 /**< ZHOU*/,
    NO_BOUNDARYFORCE = 1 /**< NO_BOUNDARYFORCE*/
  };

  /** Constructor, initializes the class  */
  BoundaryForceConfiguration() : _insertionType(NO_BOUNDARYFORCE), _density(0.0), _temperature(0.0), _isValid(true) {}

  /** Destructor */
  virtual ~BoundaryForceConfiguration() {}

  /** parseSubtag
   * 	@param node
   */
  void parseSubtag(tinyxml2::XMLElement *node) {
    std::string value;
    tarch::configuration::ParseConfiguration::readStringMandatory(value, node, "type");
    const std::string boundaries[6] = {"west", "east", "south", "north", "bottom", "top"};

    if (value == "zhou-boundary-force") {
      _insertionType = ZHOU;
    } else if (value == "none") {
      _insertionType = NO_BOUNDARYFORCE;
    } else {
      std::cout << "ERROR coupling::BoundaryForceConfiguration: Wrong insertion type!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }

    // for all boundary force types: read the sides at which boundary force
    // shall be applied
    // this may not be useful for NO_BOUNDARYFORCE; however, once read, the user
    // can still use the flags,
    // even in the NO_BOUNDARYFORCE case, to e.g. determine a setup with
    // boundary forcing using external implementations
    for (unsigned int d = 0; d < 2 * dim; d++) {
      tarch::configuration::ParseConfiguration::readBoolMandatory(_boundary[d], node, boundaries[d]);
    }

    if (_insertionType == ZHOU) {
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_density, node, "density");
      if (_density <= 0.0) {
        std::cout << "ERROR coupling::BoundaryForceConfiguration: density is "
                     "smaller than zero! Density="
                  << _density;
        _isValid = false;
        exit(EXIT_FAILURE);
      }

      tarch::configuration::ParseConfiguration::readDoubleMandatory(_temperature, node, "temperature");
      if (_temperature <= 0.0) {
        std::cout << "ERROR coupling::BoundaryForceConfiguration: temperature "
                     "is smaller than zero! Temperature="
                  << _temperature;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
    }
  }

  /** Returns name of xml tag that is associated to the configuration.
   * 	@return name of xml tag that is associated to the configuration
   */
  std::string getTag() const { return "boundary-force"; }

  /** checks if the configuration is valid. This operation usually fails, if
e.g.
         *	1. parseSubtag() hasn't been called, i.e. configuration has not
been used, or
     *  2. parseSubtag() failed due to a wrong file.
         *  3. If a tag ain't optional and parseSubtag() was not called (first
case)
         * 	@return _isValid
     */
  bool isValid() const { return _isValid; }

  /** Returns boundary force type.
   * 	@return _insertionType
   */
  const BoundaryForceType &getBoundaryForceType() const { return _insertionType; }

  /** Returns boundary ??
   * 	@return _boundary
   */
  const tarch::la::Vector<2 * dim, bool> &getBoundary() const { return _boundary; }

  /** Returns boundary force configuration.
   * 	@tparam LinkedCell type of the cell
   * 	@param mdSolverInterface
   * 	@return boundary force config
   */
  template <class LinkedCell>
  coupling::BoundaryForceController<LinkedCell, dim> *
  interpreteConfiguration(coupling::interface::MDSolverInterface<LinkedCell, dim> *const mdSolverInterface) const {
    if (_insertionType == ZHOU) {
      if (dim != 3) {
        std::cout << "ERROR "
                     "coupling::configurations::BoundaryForceConfiguration::interpret"
                     "eConfiguration(): Zhou boundary force only valid in 3D!"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      return new coupling::ZhouBoundaryForceController<LinkedCell, dim>(_density, _temperature, _boundary, mdSolverInterface);
    } else if (_insertionType == NO_BOUNDARYFORCE) {
      return new coupling::NoBoundaryForce<LinkedCell, dim>(mdSolverInterface);
    }
    return NULL;
  }

private:
  BoundaryForceType _insertionType;
  double _density;
  double _temperature;
  tarch::la::Vector<2 * dim, bool> _boundary;
  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_BOUNDARYFORCECONFIGURATION_H_
