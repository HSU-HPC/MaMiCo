#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

fileBegin = 'Molecules_0__'
fileEnd = '.vtk'

for i in range(0,340,5):
    #x = np.genfromtxt("CheckpointSimpleMD_10000_reflecting_0_0.checkpoint", delimiter=' ', usecols=(6), dtype="float", skip_header=1)
    #y = np.genfromtxt("CheckpointSimpleMD_10000_reflecting_0_0.checkpoint", delimiter=' ', usecols=(7), dtype="float", skip_header=1)
    z = np.genfromtxt(fileBegin+str(i)+fileEnd, delimiter=' ', usecols=(2), dtype="float",skip_header=6, skip_footer=2702049)
    print(np.max(z))
