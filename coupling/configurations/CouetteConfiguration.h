// Copyright (C) 2016 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _COUETTECONFIG_H_
#define _COUETTECONFIG_H_

#include "tarch/la/Vector.h"

namespace coupling {
namespace configurations {
struct CouetteConfig;
}
} // namespace coupling
/** Configuration parameters for Couette flow scenario
 *	@brief Configuration parameters for Couette flow scenario
 *
 *
 *	@author Philipp Neumann
 */
struct coupling::configurations::CouetteConfig {

private:
  // CouetteConfig(){}

public:
  /** Macroscopic solver types that are supported.
   *	@enum MacroSolverType
   */
  enum MacroSolverType {
    COUETTE_ANALYTICAL = 0 /**< Analytic solver @note only supports flow in
x-direction, is only active on rank 0, and that it does not model viscous shear
forces in the oscillating wall case.*/
    ,
    COUETTE_LB = 1 /**< Lattice Boltzmann solver @note can support flow in all
     direction (3D) and can be parallelized*/
    ,
    COUETTE_FD = 2 /**< Finite difference solver*/
  };
  /** Microscopic solver types that are supported.
   *	@enum MicroSolverType
   */
  enum MicroSolverType {
    SIMPLEMD = 0 /**< simple MD*/,
    SYNTHETIC = 1 /**< data from macroscopic plus gaussian noise*/
  };

  /** @brief operator overloading; declare a configuration (this)
   * 	@param other
   */
  CouetteConfig &operator=(const CouetteConfig &other) {

    this->channelheight = other.channelheight;
    this->wallInitCycles = other.wallInitCycles;
    this->wallVelocity = other.wallVelocity;
    this->wallOscillations = other.wallOscillations;
    this->couplingCycles = other.couplingCycles;
    this->md2Macro = other.md2Macro;
    this->macro2Md = other.macro2Md;
    this->computeSNR = other.computeSNR;
    this->twoWayCoupling = other.twoWayCoupling;
    this->filterInitCycles = other.filterInitCycles;
    this->csvEveryTimestep = other.csvEveryTimestep;
    this->maSolverType = other.maSolverType;
    this->miSolverType = other.miSolverType;
    this->density = other.density;
    this->kinVisc = other.kinVisc;
    this->lbNumberProcesses = other.lbNumberProcesses;
    this->plotEveryTimestep = other.plotEveryTimestep;
    this->initAdvanceCycles = other.initAdvanceCycles;
    this->temp = other.temp;
    this->equSteps = other.equSteps;
    this->totalNumberMDSimulations = other.totalNumberMDSimulations;
    this->noiseSigma = other.noiseSigma;
    this->twsLoop = other.twsLoop;
    this->twsLoopMax = other.twsLoopMax;
    this->twsLoopMin = other.twsLoopMin;
    this->twsLoopStep = other.twsLoopStep;

    return *this;
  }

  /** @brief creates CouetteConfig if all elements exist and can be read
   * 	@param filename
   */
  static CouetteConfig parseCouetteConfiguration(const std::string &filename) {
    CouetteConfig _cfg;

    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement *node = NULL;
    conffile.LoadFile("couette.xml");
    node = conffile.FirstChildElement("couette-test");
    if (node == NULL) {
      std::cout << "Could not read input file couette.xml: missing element "
                   "<couette-test>"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement *n_mamico = node->NextSiblingElement();
    if (n_mamico == NULL) {
      std::cout << "Could not read input file couette.xml: missing element <mamico>" << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n_md = n_mamico->NextSiblingElement();
    if (n_md == NULL) {
      std::cout << "Could not read input file couette.xml: missing element "
                   "<molecular-dynamics>"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n_fp = n_md->NextSiblingElement();
    if (n_fp == NULL) {
      std::cout << "Could not read input file couette.xml: missing element "
                   "<filter-pipeline>"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tinyxml2::XMLElement *n_unexpected = n_fp->NextSiblingElement();
    if (n_unexpected == NULL) {
      std::cout << "Could not read input file couette.xml: unknown element " << n_unexpected->Name() << std::endl;
      exit(EXIT_FAILURE);
    }

    tinyxml2::XMLElement *subtag = node->FirstChildElement("domain");
    if (subtag == NULL) {
      std::cout << "Could not read input file couette.xml: Missing subtag: domain" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.channelheight, subtag, "channelheight");
    tarch::configuration::ParseConfiguration::readVector<3, double>(_cfg.wallVelocity, subtag, "wall-velocity");
    _cfg.wallInitCycles = 0;
    tarch::configuration::ParseConfiguration::readIntOptional(_cfg.wallInitCycles, subtag, "wall-init-cycles");
    if (_cfg.wallInitCycles > 0)
      tarch::configuration::ParseConfiguration::readVector<3, double>(_cfg.wallInitVelocity, subtag, "wall-init-velocity");
    _cfg.wallOscillations = 0;
    tarch::configuration::ParseConfiguration::readDoubleOptional(_cfg.wallOscillations, subtag, "wall-oscillations");

    subtag = node->FirstChildElement("coupling");
    if (subtag == NULL) {
      std::cout << "Could not read input file couette.xml: Missing subtag: coupling" << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.couplingCycles, subtag, "coupling-cycles");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.twoWayCoupling, subtag, "two-way-coupling");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.md2Macro, subtag, "send-from-md-to-macro");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.macro2Md, subtag, "send-from-macro-to-md");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.filterInitCycles, subtag, "filter-init-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep, subtag, "write-csv-every-timestep");
    tarch::configuration::ParseConfiguration::readBoolMandatory(_cfg.computeSNR, subtag, "compute-snr");

    subtag = node->FirstChildElement("microscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file couette.xml: Missing subtag: "
                   "microscopic-solver"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    if (type == "md") {
      _cfg.miSolverType = SIMPLEMD;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp, subtag, "temperature");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps, subtag, "equilibration-steps");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
      if (_cfg.totalNumberMDSimulations < 1) {
        std::cout << "Could not read input file couette.xml: "
                     "number-md-simulations < 1"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if (type == "synthetic") {
      _cfg.miSolverType = SYNTHETIC;
      tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.noiseSigma, subtag, "noise-sigma");
      _cfg.totalNumberMDSimulations = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
    } else {
      std::cout << "Could not read input file couette.xml: Unknown microscopic "
                   "solver type!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.density, subtag, "density");

    subtag = node->FirstChildElement("macroscopic-solver");
    if (subtag == NULL) {
      std::cout << "Could not read input file couette.xml: Missing subtag: "
                   "macroscopic-solver"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    _cfg.lbNumberProcesses = tarch::la::Vector<3, unsigned int>(1);
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    if (type == "lb") {
      _cfg.maSolverType = COUETTE_LB;
      tarch::configuration::ParseConfiguration::readVector<3, unsigned int>(_cfg.lbNumberProcesses, subtag, "number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.plotEveryTimestep, subtag, "plot-every-timestep");
    } else if (type == "fd") {
      _cfg.maSolverType = COUETTE_FD;
      tarch::configuration::ParseConfiguration::readVector<3, unsigned int>(_cfg.lbNumberProcesses, subtag, "number-of-processes");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.plotEveryTimestep, subtag, "plot-every-timestep");
    } else if (type == "analytical") {
      _cfg.maSolverType = COUETTE_ANALYTICAL;
      if (!(_cfg.wallVelocity[1] == 0.0 && _cfg.wallVelocity[2] == 0.0)) {
        std::cout << "analytic solver only supports flow in x-direction" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cout << "Could not read input file couette.xml: Unknown macroscopic "
                   "solver type!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    double vis;
    tarch::configuration::ParseConfiguration::readDoubleMandatory(vis, subtag, "viscosity");
    _cfg.kinVisc = vis / _cfg.density;
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.initAdvanceCycles, subtag, "init-advance-cycles");

    subtag = node->FirstChildElement("tws-loop");
    if (subtag == NULL)
      _cfg.twsLoop = false;
    else {
      _cfg.twsLoop = true;
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMin, subtag, "min");
      tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.twsLoopMax, subtag, "max");
      _cfg.twsLoopStep = 1;
      tarch::configuration::ParseConfiguration::readIntOptional(_cfg.twsLoopStep, subtag, "step");
    }

    if (_cfg.miSolverType == SYNTHETIC) {
      if (_cfg.md2Macro || _cfg.macro2Md || _cfg.totalNumberMDSimulations > 1 || _cfg.lbNumberProcesses[0] != 1 || _cfg.lbNumberProcesses[1] != 1 ||
          _cfg.lbNumberProcesses[2] != 1) {
        std::cout << "Invalid configuration: Synthetic MD runs sequentially on "
                     "rank 0 only. "
                  << "It does neither support parallel communication nor "
                     "multi-instance sampling"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (_cfg.maSolverType == COUETTE_ANALYTICAL) {
      if (_cfg.twoWayCoupling) {
        std::cout << "Invalid configuration: COUETTE_ANALYTICAL does not "
                     "support twoWayCoupling"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    return _cfg;
  }

  /** channel is always expected to have origin at (0.0,0.0,0.0) and to be cubic
   * (MD 30: 50.0, MD 60: 100.0, MD 120: 200.0) */
  double channelheight;
  /** velocity of moving wall (lower boundary moves)*/
  tarch::la::Vector<3, double> wallVelocity;
  int wallInitCycles;
  tarch::la::Vector<3, double> wallInitVelocity;
  double wallOscillations;
  /** number of coupling cycles, that is continuum time steps; MD/DPD: 1000*/
  int couplingCycles;
  bool md2Macro, macro2Md, computeSNR, twoWayCoupling;
  int filterInitCycles;
  int csvEveryTimestep;
  MacroSolverType maSolverType;
  MicroSolverType miSolverType;
  double density, kinVisc;
  /** only for LB couette solver: number of processes */
  tarch::la::Vector<3, unsigned int> lbNumberProcesses;
  /** only for LB couette solver: VTK plotting per time step*/
  int plotEveryTimestep;
  int initAdvanceCycles;
  double temp;
  int equSteps;
  int totalNumberMDSimulations;
  double noiseSigma;
  bool twsLoop;
  int twsLoopMin, twsLoopMax, twsLoopStep;
};

#endif