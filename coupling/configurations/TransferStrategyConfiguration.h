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
template <unsigned int dim> class TransferStrategyConfiguration;
}
}

/** transfer strategy configuration, i.e. algorithm/combin. of quantity transfer
 * steps and quantity interpretation (e.g. momentum vs. velocity). Derive from
 * the class tarch::configuration::Configuration
 *	@brief transfer strategy configuration, i.e. algorithm/combin. of quantity
 * transfer steps and quantity interpretation (e.g. momentum vs. velocity).
 *	@tparam dim Number of dimensions; it can be 1, 2 or 3
 *  @author Philipp Neumann
 */
template <unsigned int dim>
class coupling::configurations::TransferStrategyConfiguration
    : public tarch::configuration::Configuration {
public:
  /** transfer strategy types that are implemented.
 	 *	@enum StrategyType
 	 */
  enum StrategyType {
    DirectTransferStrategy = 0 /**< DirectTransferStrategy*/,
    DifferenceTransferStrategy = 1 /**< DifferenceTransferStrategy*/,
    TransferStrategy4FluxCoupling = 2 /**< TransferStrategy4FluxCoupling*/,
    TransferStrategy4SchwarzCoupling =
        3 /**< TransferStrategy4SchwarzCoupling*/,
    AveragingTransferStrategy = 4 /**< AveragingTransferStrategy*/,
    TransferStrategy4NieCoupling = 5 /**< TransferStrategy4NieCoupling*/
  };

  /** Constructor, initializes the class  */
  TransferStrategyConfiguration()
      : _type(DirectTransferStrategy), _massFluxBoundary(false),
        _isValid(true) {}

  /** Destructor */
  virtual ~TransferStrategyConfiguration() {}

  /** parseSubtag
	 * 	@param node
     */
  void parseSubtag(tinyxml2::XMLElement *node) {
    std::string value;
    tarch::configuration::ParseConfiguration::readStringMandatory(value, node,
                                                                  "type");
    if (value == "direct-transfer") {
      _type = DirectTransferStrategy;
    } else if (value == "difference-transfer") {
      _type = DifferenceTransferStrategy;
    } else if (value == "schwarz-transfer") {
      _type = TransferStrategy4SchwarzCoupling;
    } else if (value == "nie-transfer") {
      _type = TransferStrategy4NieCoupling;
    } else if (value == "averaging") {
      _type = AveragingTransferStrategy;
    } else {
      std::cout << "ERROR coupling::TransferStrategyConfiguration: Wrong "
                   "insertion type!" << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }

    if (_type == TransferStrategy4NieCoupling) {
      const std::string boundaries[6] = { "mass-flux-west", "mass-flux-east",
                                          "mass-flux-south", "mass-flux-north",
                                          "mass-flux-bottom", "mass-flux-top" };
      for (unsigned int d = 0; d < 2 * dim; d++) {
        tarch::configuration::ParseConfiguration::readBoolMandatory(
            _massFluxBoundary[d], node, boundaries[d]);
      }
      // by default: no shifting of time interval
      _shiftTimestep = 0.0;
      tarch::configuration::ParseConfiguration::readDoubleOptional(
          _shiftTimestep, node, "shift-by-timesteps");
      if (_shiftTimestep < 0.0 || _shiftTimestep > 1.0) {
        std::cout << "Warning "
                     "coupling::configurations::TransferStrategyConfiguration: "
                     "shift-by-timesteps=" << _shiftTimestep
                  << "; typical values range between 0 and 1!" << std::endl;
      }
    }
  }

  /** Returns name of xml tag that is associated to the configuration.
	 * 	@return name of xml tag that is associated to the configuration
     */
  std::string getTag() const { return "transfer-strategy"; }

  /** checks if the configuration is valid. This operation usually fails, if
e.g.
	 *	1. parseSubtag() hasn't been called, i.e. configuration has not been used,
or
     *  2. parseSubtag() failed due to a wrong file.
	 *  3. If a tag ain't optional and parseSubtag() was not called (first case)
	 * 	@return _isValid
     */
  bool isValid() const { return _isValid; }

  /** Returns transfer strategy configuration.
	 * 	@tparam LinkedCell type of the cell
	 * 	@param mdSolverInterface
	 * 	@param indexConversion
	 * 	@param numberOfMDTimesteps
	 * 	@return transfer strategy config
     */
  template <class LinkedCell>
  coupling::transferstrategies::TransferStrategy<LinkedCell, dim> *
  interpreteConfiguration(coupling::interface::MDSolverInterface<
                              LinkedCell, dim> *const mdSolverInterface,
                          const coupling::IndexConversion<dim> &indexConversion,
                          unsigned int numberOfMDTimesteps) const {
    if (_type == DirectTransferStrategy) {
      return new coupling::transferstrategies::DirectTransferStrategy<
          LinkedCell, dim>(mdSolverInterface, indexConversion);
    } else if (_type == DifferenceTransferStrategy) {
      return new coupling::transferstrategies::DifferenceTransferStrategy<
          LinkedCell, dim>(mdSolverInterface, indexConversion,
                           numberOfMDTimesteps);
    } else if (_type == TransferStrategy4SchwarzCoupling) {
      return new coupling::transferstrategies::TransferStrategy4SchwarzCoupling<
          LinkedCell, dim>(mdSolverInterface, indexConversion,
                           numberOfMDTimesteps);
    } else if (_type == TransferStrategy4NieCoupling) {
      return new coupling::transferstrategies::TransferStrategy4NieCoupling<
          LinkedCell, dim>(mdSolverInterface, indexConversion,
                           numberOfMDTimesteps, _shiftTimestep,
                           _massFluxBoundary);
    } else if (_type == AveragingTransferStrategy) {
      return new coupling::transferstrategies::AveragingTransferStrategy<
          LinkedCell, dim>(mdSolverInterface, indexConversion);
    } else {
      return NULL;
    }
  }

  /** Returns the transfer strategy type.
	 * 	@return _type
     */
  StrategyType getStrategyType() const { return _type; }

private:
  StrategyType _type;
  tarch::la::Vector<2 * dim, bool>
      _massFluxBoundary; // true in each component, if one of the 2*dim
                         // boundaries allows for mass flux
  double _shiftTimestep; // used for Nie coupling: time interval by which the
                         // evaluation of the continuum flow field should be
                         // shifted. See also TransferStrategy4NieCoupling for
                         // more details.

  bool _isValid;
};

#endif // _MOLECULARDYNAMICS_COUPLING_CONFIGURATIONS_STRATEGYCONFIGURATION_H_
