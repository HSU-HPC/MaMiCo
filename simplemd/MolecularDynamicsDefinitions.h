// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULARDYNAMICSDEFINITIONS_H_
#define _MOLECULARDYNAMICS_MOLECULARDYNAMICSDEFINITIONS_H_

#include "simplemd/MolecularDynamicsUserInput.h"
#include "tarch/la/Vector.h"

#if (MD_DIM == 1)
#define MD_LINKED_CELL_NEIGHBOURS 2
#elif (MD_DIM == 2)
#define MD_LINKED_CELL_NEIGHBOURS 8
#elif (MD_DIM == 3)
#define MD_LINKED_CELL_NEIGHBOURS 26
#else
#error "Only 1D/2D/3D supported!"
#endif

#define MD_PI 3.1415926535897932384626433832795028841972
#define MD_TOLERANCE 1e-13

namespace simplemd {
// static const unsigned int NUMBER_MOLECULES_PER_BLOCK_ALLOCATION = 1;

// PERIODIC_BOUNDARY: Ghost linked cells are filled by respective molecular
// information from the cells on opposite end of box
// OPEN_BOUNDARY: Additional force term is applied to molecules near this
// boundary to exert correct hydrodynamic pressure
//                -> currently not supported!
// PARALLEL_BOUNDARY: Similar to PERIODIC_BOUNDARY, but molecules need to be
// fetched from neighboring process
// GEOMETRY_BOUNDARY: No molecules are assumed to be in these cells (geometry
// needs to be included via the respective
//                    geometrical forms). -> currently not supported!
enum BoundaryType { NO_BOUNDARY = 0, PERIODIC_BOUNDARY = 1, OPEN_BOUNDARY = 2, PARALLEL_BOUNDARY = 3, GEOMETRY_BOUNDARY = 4, REFLECTING_BOUNDARY = 5 };

// static const tarch::la::Vector<MD_DIM,double> zero(0.0);

} // namespace simplemd

#endif // _MOLECULARDYNAMICS_MOLECULARDYNAMICSDEFINITIONS_H_
