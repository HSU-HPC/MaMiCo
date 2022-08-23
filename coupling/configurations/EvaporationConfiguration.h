#ifndef _EVAPORATIONCONFIG_H_
#define _EVAPORATIONCONFIG_H_

#include "tarch/la/Vector.h"

namespace coupling {
namespace configurations {
struct EvaporationConfig;
}
}

struct coupling::configurations::EvaporationConfig {

public:
  static EvaporationConfig parseConfiguration(const std::string& filename) {
    EvaporationConfig _cfg;

    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement* node = NULL;
    conffile.LoadFile("evaporation.xml");
    node = conffile.FirstChildElement("scenario");
    tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.couplingCycles, subtag, "coupling-cycles");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep, subtag, "write-csv-every-timestep");
    subtag = node->FirstChildElement("microscopic-solver");
    std::string type;
    tarch::configuration::ParseConfiguration::readStringMandatory(type, subtag, "type");
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp, subtag, "temperature");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps, subtag, "equilibration-steps");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
    return _cfg;
  }

  int couplingCycles;
  int csvEveryTimestep;
  int equSteps;
  double temp;
  int totalNumberMDSimulations;
};

#endif
