// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_TARCHDEFINITIONS_H_
#define _TARCH_TARCHDEFINITIONS_H_

#define TARCH_YES 1
#define TARCH_NO 0

#ifdef TarchParallel
#define TARCH_PARALLEL TARCH_YES
#else
#define TARCH_PARALLEL TARCH_NO
#endif

#ifdef TarchDebug
#define TARCH_DEBUG TARCH_YES
#else
#define TARCH_DEBUG TARCH_NO
#endif

#endif // _TARCH_TARCHDEFINITIONS_H_
