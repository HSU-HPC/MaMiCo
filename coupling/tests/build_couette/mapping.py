#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t_start = 0
t_end = 34
t_step = 2
rho_mean = np.zeros(np.int(t_end/t_step))
c = 0
for t in np.arange(t_start, t_end, t_step):
    rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(t)+".csv", delimiter=';', usecols=(6), dtype="float")
    rho_mean[c] = np.mean(rho/15.625)
    c=c+1
fig = plt.figure()
plt.plot(rho_mean[:c])
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
