#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

cellSize = 200*200*2.5

fig = plt.figure()

for i in range(1,11,1):
    file = 'CouetteAvgMultiMDCells_0_'+str(i)+'.csv'
    u = np.genfromtxt(file, usecols=(3), delimiter=';', dtype="double")
    v = np.genfromtxt(file, usecols=(4), delimiter=';', dtype="double")
    w = np.genfromtxt(file, usecols=(5), delimiter=';', dtype="double")
    rho = np.genfromtxt(file, usecols=(6), delimiter=';', dtype="double")
    print(np.sum(rho))
    rho /= cellSize

    #plt.plot(np.arange(169)*2.5, u, label='u '+str(i))
    #plt.plot(np.arange(169)*2.5, v, label='v '+str(i))
    #plt.plot(np.arange(169)*2.5, w, label='w '+str(i))
    plt.plot(np.arange(169)*2.5, rho, label='rho '+str(i))
plt.xlabel("z")
plt.legend()
plt.xlim(0,422.5)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
