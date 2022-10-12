// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace tarch {
namespace logging {
class Logger {
public:
  Logger() : Logger("") {}

  Logger(std::string name) {
    _spdlogger = spdlog::get(name);
    if (_spdlogger == NULL) _spdlogger = spdlog::stdout_color_mt(name);
    _spdlogger->set_level(spdlog::level::debug);
    _spdlogger->set_pattern(DEFAULT_PATTERN);
  }

  template <typename... Args> 
  void info(spdlog::format_string_t<Args...> fmt, Args &&... args) const {
     _spdlogger->info(fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> 
  void debug(spdlog::format_string_t<Args...> fmt, Args &&... args) const {
    _spdlogger->debug(fmt, std::forward<Args>(args)...);
  }

  template <typename... Args> 
  void error(spdlog::format_string_t<Args...> fmt, Args &&... args) const {
    _spdlogger->error(fmt, std::forward<Args>(args)...);
  }

private:
  static const std::string DEFAULT_PATTERN;
  std::shared_ptr<spdlog::logger> _spdlogger;
};
} // namespace tarcg
} // namespace logging

const std::string tarch::logging::Logger::DEFAULT_PATTERN = "[%H:%M:%S %^%l%$] MaMiCo %n: %v";