#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

cellSize = 200*200*.25

fig = plt.figure()

for i in range(49,250,50):
    file = 'EvaporationZAvgBin_0_'+str(i)+'.csv'
    u = np.genfromtxt(file, usecols=(1), delimiter=';', dtype="double")
    v = np.genfromtxt(file, usecols=(2), delimiter=';', dtype="double")
    w = np.genfromtxt(file, usecols=(3), delimiter=';', dtype="double")
    rho = np.genfromtxt(file, usecols=(4), delimiter=';', dtype="double")/cellSize
    T = np.genfromtxt(file, usecols=(5), delimiter=';', dtype="double")
    wPlus = np.genfromtxt(file, usecols=(8), delimiter=';', dtype="double")
    wMinus = np.genfromtxt(file, usecols=(12), delimiter=';', dtype="double")
    rhoPlus = np.genfromtxt(file, usecols=(9), delimiter=';', dtype="double")/cellSize
    rhoMinus = np.genfromtxt(file, usecols=(13), delimiter=';', dtype="double")/cellSize

    #plt.plot(np.arange(1680)*.25, u, label='u '+str(i))
    #plt.plot(np.arange(1680)*.25, v, label='v '+str(i))
    #plt.plot(np.arange(1680)*.25, w, label='w '+str(i))
    #plt.plot(np.arange(1680)*.25, rho, label='rho '+str(i))
    plt.plot(np.arange(1680)*.25, T, label='T '+str(i))
    #plt.plot(np.arange(1680)*.25, rho*w, label='j_p '+str(i))
    #plt.plot(np.arange(1680)*.25, rhoPlus, label='rhoPlus '+str(i))
    #plt.plot(np.arange(1680)*.25, rhoMinus, label='rhoMins '+str(i))
    #plt.plot(np.arange(1680)*.25, rhoPlus+rhoMinus, label='rho+ '+str(i))
    #plt.plot(np.arange(1680)*.25, wPlus, label='wPlus '+str(i))
    #plt.plot(np.arange(1680)*.25, wMinus, label='wMinus '+str(i))
    #plt.plot(np.arange(1680)*.25, wPlus*rhoPlus, label='jPlus '+str(i))
    #plt.plot(np.arange(1680)*.25, -wMinus*rhoMinus, label='jMinus '+str(i))
    #plt.plot(np.arange(1680)*.25, wPlus*rhoPlus+wMinus*rhoMinus, label='jMinus '+str(i))
plt.xlabel("z")
plt.legend()
plt.xlim(0,422.5)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
