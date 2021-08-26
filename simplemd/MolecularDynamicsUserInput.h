// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _MOLECULARDYNAMICS_MOLECULARDYNAMICSUSERINPUT_H_
#define _MOLECULARDYNAMICS_MOLECULARDYNAMICSUSERINPUT_H_

#define MD_YES 1
#define MD_NO 0

// USER INPUT --------------------------
// Dimension of simulation (1D=1,2D=2, 3D=3)
#if defined(MDDim2)
#define MD_DIM 2
#elif defined(MDDim3)
#define MD_DIM 3
#else
#error "MDDim2 or MDDim3"
#endif

// MPI-Parallelisation (MD_YES/ MD_NO)
#ifdef MDParallel
#define MD_PARALLEL MD_YES
#else
#define MD_PARALLEL MD_NO
#endif

// debug mode (MD_YES/MD_NO). Prints very lot of stuff to the screen. Should only be used
// during debugging.
#ifdef MDDebug
#define MD_DEBUG MD_YES
#else
#define MD_DEBUG MD_NO
#endif

#ifdef MDError
#define MD_ERROR MD_YES
#else
// error mode (MD_YES/MD_NO). Performs error checks during the computation.
#define MD_ERROR MD_NO
#endif

// OpenMP-Parallelisation (MD_YES/MD_NO)
#define MD_OPENMP MD_NO

// TEST_TCHIPEV (MD_YES/MD_NO). If no - two communication stages, if yes - one communication stage.
#define PARALLEL_GDB_DEBUG MD_NO
#endif // _MOLECULARDYNAMICS_MOLECULARDYNAMICSUSERINPUT_H_

