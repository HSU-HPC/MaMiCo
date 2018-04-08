// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_COUPLING_TESTS_TEST_H_
#define _MOLECULARDYNAMICS_COUPLING_TESTS_TEST_H_

#include <string>
#include <iostream>

class Test {
  public:
    Test(std::string testname): _testname(testname){
      std::cout << "Run " << testname << "..." << std::endl;
    }
    virtual ~Test(){
      std::cout << "Shut down " << _testname << std::endl;
    }

    virtual void run() = 0;

  private:
    const std::string _testname;
};

#endif // _MOLECULARDYNAMICS_COUPLING_TESTS_TEST_H_
