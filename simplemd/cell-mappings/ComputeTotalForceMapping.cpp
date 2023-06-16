// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/ComputeTotalForceMapping.h"

void simplemd::cellmappings::ComputeTotalForceMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  // iterate over all molecules from this cell
  for (std::list<Molecule*>::iterator it = cell.begin(); it != cell.end(); it++) {
    // get the molecule pointer
    Molecule* molecule = (*it);

    // calculate the total force from the partial forces
    double position = molecule->getPosition()[_dimension - 1];
    double weight = calculateWeight(position);
    tarch::la::Vector<MD_DIM, double> totalForce = weight * molecule->getConstThreeBodyForce() + (1.0 - weight) * molecule->getConstTwoBodyForce();
    molecule->setForce(totalForce);
    molecule->setTwoBodyForce(_zero);
    molecule->setThreeBodyForce(_zero);
  }
}
