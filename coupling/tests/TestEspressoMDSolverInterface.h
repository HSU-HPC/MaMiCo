// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "global.hpp"

class TestFactory {
public:
  static void getReferenceSetup() {
    box_l[0] = box_l[1] = box_l[2] = 10.0;
    dd.cell_grid[0] = dd.cell_grid[1] = dd.cell_grid[2] = 4;
    double number_density = 0.6;
    n_part = number_density * box_l[0] * box_l[1] * box_l[2];
  }
} class TestCase {
public:
  void run() {
    TestFactory _factory.getReferenceSetup();
    //TODO rest of test case
  }
}
