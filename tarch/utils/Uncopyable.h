// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_UTILS_UNCOPYABLE_H_
#define _TARCH_UTILS_UNCOPYABLE_H_

namespace tarch {
namespace utils { class Uncopyable; }
}

// class to prevent copying; just inherit "privately" from this and copying will
// not work anymore (see Item6 of Effective C++, S. Meyers)
class tarch::utils::Uncopyable {
protected:
  Uncopyable() {}
  ~Uncopyable() {}

private:
  Uncopyable(const tarch::utils::Uncopyable &);
  Uncopyable &operator=(const tarch::utils::Uncopyable &);
};

#endif //
