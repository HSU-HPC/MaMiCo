#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

cellSize = 200*200*.25

fig = plt.figure()

for i in range(2999,3500,50):
    file = 'EvaporationZAvgBin_'+str(i)+'.csv'
    u = np.genfromtxt(file, usecols=(1), delimiter=';', dtype="double")
    v = np.genfromtxt(file, usecols=(2), delimiter=';', dtype="double")
    w = np.genfromtxt(file, usecols=(3), delimiter=';', dtype="double")
    rho = np.genfromtxt(file, usecols=(4), delimiter=';', dtype="double")/cellSize
    T = np.genfromtxt(file, usecols=(5), delimiter=';', dtype="double")

    #plt.plot(np.arange(1680)*.25, u, label='u '+str(i))
    #plt.plot(np.arange(1680)*.25, v, label='v '+str(i))
    #plt.plot(np.arange(1680)*.25, w, label='w '+str(i))
    plt.plot(np.arange(1680)*.25, rho, label='rho '+str(i))
    #plt.plot(np.arange(1680)*.25, T, label='T '+str(i))
    #plt.plot(np.arange(1680)*.25, rho*w, label='j_p '+str(i))
plt.xlabel("z")
plt.legend()
plt.xlim(0,422.5)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
