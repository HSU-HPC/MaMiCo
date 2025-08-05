// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#include "simplemd/cell-mappings/DeleteMoleculesMapping.h"

void simplemd::cellmappings::DeleteMoleculesMapping::handleCell(LinkedCell& cell, const unsigned int& cellIndex) {
  // clear cell list
  cell.clear();
}
