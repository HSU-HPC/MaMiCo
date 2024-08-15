// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder

#include "DomainDecompConfiguration.h"
#include "tarch/configuration/ParseConfiguration.h"

const std::string simplemd::configurations::DomainDecompConfiguration::DECOMP_TYPE("decomposition-type");
const std::string simplemd::configurations::DomainDecompConfiguration::DEFAULT_DECOMP("default");
const std::string simplemd::configurations::DomainDecompConfiguration::STATIC_IRREG_RECT_GRID("static-irreg-rect-grid");
const std::string simplemd::configurations::DomainDecompConfiguration::AXES[MD_DIM]
#if MD_DIM == 1
    = {"x"};
#elif MD_DIM == 2
    = {"x", "y"};
#elif MD_DIM == 3
    = {"x", "y", "z"};
#endif

simplemd::configurations::DomainDecompConfiguration::DomainDecompConfiguration() {
  _decompType = DecompositionType::DEFAULT;
  _isDefined = false;
  _isValid = true;
  for (int i = 0; i < MD_DIM; i++) {
    _subdomainWeights[i] = {};
  }
}

void simplemd::configurations::DomainDecompConfiguration::parseSubtag(tinyxml2::XMLElement* node) {
  _isDefined = true;
  std::string stringBuf;

  tarch::configuration::ParseConfiguration::readStringMandatory(stringBuf, node, DECOMP_TYPE);
  if (stringBuf == DEFAULT_DECOMP)
    _decompType = DecompositionType::DEFAULT;
  else if (stringBuf == STATIC_IRREG_RECT_GRID)
    _decompType = DecompositionType::STATIC_IRREG_RECT_GRID;
  else {
    std::cout << "ERROR: given decomposition-type not supported!" << std::endl;
    _isValid = false;
    exit(EXIT_FAILURE);
  }

  switch (_decompType) {
  case DecompositionType::DEFAULT:
    /* nop */
    break;

  case DecompositionType::STATIC_IRREG_RECT_GRID: {
    std::string weightsBuf;

    for (unsigned int d = 0; d < MD_DIM; d++) {
      tarch::configuration::ParseConfiguration::readStringMandatory(weightsBuf, node, AXES[d]);
      if (weightsBuf.length() == 0) {
        std::cout << "ERROR: empty weights for axis " << AXES[d] << std::endl;
        _isValid = false;
        exit(EXIT_FAILURE);
      }
      _subdomainWeights[d] = getWeightsFromString(weightsBuf);
    }
  } break;

  default:
    break;
  }
}

std::vector<unsigned int> simplemd::configurations::DomainDecompConfiguration::getWeightsFromString(std::string weights) {
  std::vector<unsigned int> result;
  unsigned int temp;
  std::stringstream ss(weights);
  // Taken partially from ls1-mardyn::StaticIrregDomainDecomposition.cpp
  // Parse the weights, until the stringstream has chars and extraction
  // doesn't fail, and no EOF or linebreaks etc
  while (ss.good()) {
    while (ss.peek() == ';' || ss.peek() == ' ') // skip semicolons and spaces
      ss.ignore();
    // Extraction from stream into int type fails if token is not an int
    // We check for this failure, and additionally check for positive
    // integer
    if (!(ss >> temp) || temp <= 0) {
      std::cout << "Weights (" << weights
                << ") have a non-natural number! Only integer weights > "
                   "0 allowed, please check XML file!"
                << std::endl;
      _isValid = false;
      exit(EXIT_FAILURE);
    }
    result.push_back(temp);
  }
  return result;
}

std::string simplemd::configurations::DomainDecompConfiguration::getTag() const { return "domain-decomp-configuration"; }

bool simplemd::configurations::DomainDecompConfiguration::isValid() const { return _isValid; }
