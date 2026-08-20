// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

struct CouetteCommandLineOptions {
  std::string configurationFile{"couette.xml"};
  bool configurationFileSpecified{false};
  bool showHelp{false};
  std::vector<std::string> remainingArguments;
};

/**
 * Extract Couette-specific command line options while retaining options for
 * MPI and Kokkos.  The returned remaining arguments always include argv[0].
 */
inline CouetteCommandLineOptions parseCouetteCommandLine(const std::vector<std::string>& arguments) {
  CouetteCommandLineOptions options;
  if (arguments.empty())
    return options;

  options.remainingArguments.push_back(arguments.front());
  for (std::size_t index = 1; index < arguments.size(); ++index) {
    const std::string& argument = arguments[index];
    if (argument == "--help" || argument == "-h") {
      options.showHelp = true;
    } else if (argument == "--config") {
      if (++index == arguments.size() || arguments[index].empty())
        throw std::invalid_argument("--config requires a non-empty path");
      if (options.configurationFileSpecified)
        throw std::invalid_argument("--config may only be specified once");
      options.configurationFile = arguments[index];
      options.configurationFileSpecified = true;
    } else if (argument.rfind("--config=", 0) == 0) {
      const std::string configurationFile = argument.substr(std::string("--config=").size());
      if (configurationFile.empty())
        throw std::invalid_argument("--config requires a non-empty path");
      if (options.configurationFileSpecified)
        throw std::invalid_argument("--config may only be specified once");
      options.configurationFile = configurationFile;
      options.configurationFileSpecified = true;
    } else {
      options.remainingArguments.push_back(argument);
    }
  }
  return options;
}
