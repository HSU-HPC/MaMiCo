#ifndef _EVAPORATIONCONFIG_H_
#define _EVAPORATIONCONFIG_H_

#include "tarch/la/Vector.h"

namespace coupling {
namespace configurations {
struct ScenarioConfig;
}
}

struct coupling::configurations::ScenarioConfig {

public:
  static ScenarioConfig parseConfiguration(const std::string& filename) {
    ScenarioConfig _cfg;
    tinyxml2::XMLDocument conffile;
    tinyxml2::XMLElement* node = NULL;
    conffile.LoadFile(filename.c_str());
    node = conffile.FirstChildElement("scenario");
    tinyxml2::XMLElement* subtag = node->FirstChildElement("coupling");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.csvEveryTimestep, subtag, "write-csv-every-timestep");
    subtag = node->FirstChildElement("microscopic-solver");
    tarch::configuration::ParseConfiguration::readDoubleMandatory(_cfg.temp, subtag, "temperature");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.equSteps, subtag, "equilibration-steps");
    tarch::configuration::ParseConfiguration::readIntMandatory(_cfg.totalNumberMDSimulations, subtag, "number-md-simulations");
    return _cfg;
  }

  int csvEveryTimestep;
  int equSteps;
  double temp;
  int totalNumberMDSimulations;
};

#endif
