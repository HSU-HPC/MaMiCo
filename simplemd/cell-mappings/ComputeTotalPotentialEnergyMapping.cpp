// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/ComputeTotalPotentialEnergyMapping.h"

void simplemd::cellmappings::ComputeTotalPotentialEnergyMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  // iterate over all molecules from this cell
  for (std::list<Molecule*>::iterator it = cell.begin(); it != cell.end(); it++) {
    // get the molecule pointer
    Molecule* molecule = (*it);

    // calculate the total potential from the partial potentials
    double position = molecule->getPosition()[_dimension - 1];
    double weight = calculateWeight(position);
    double totalPotentialEnergy = weight * molecule->getConstThreeBodyPotentialEnergy() + (1.0 - weight) * molecule->getConstTwoBodyPotentialEnergy();
    molecule->setPotentialEnergy(totalPotentialEnergy);
    molecule->setTwoBodyPotentialEnergy(0.0);
    molecule->setThreeBodyPotentialEnergy(0.0);
  }
}