#!/usr/bin/env python3

import mamico.tarch.utils

rank = mamico.tarch.utils.initMPI()

print("rank = " + str(rank))

mamico.tarch.utils.finalizeMPI()
