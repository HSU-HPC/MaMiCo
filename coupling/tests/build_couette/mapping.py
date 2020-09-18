#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t_start = 0
t_end = 86
t_step = 2
volume = 5*5*5
rho_mean = np.zeros(np.int(t_end/t_step))
rho_left = np.zeros(np.int(t_end/t_step))
rho_middle = np.zeros(np.int(t_end/t_step))
rho_right = np.zeros(np.int(t_end/t_step))
c = 0
for t in np.arange(t_start, t_end, t_step):
    rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(t)+".csv", delimiter=';', usecols=(6), dtype="float")
    rho_mean[c] = np.mean(rho/volume)
    rho_left[c] = np.mean(rho[np.arange(732,847,12)])/volume
    rho_middle[c] = np.mean(rho[np.arange(738,847,12)])/volume
    rho_right[c] = np.mean(rho[np.arange(743,854,12)])/volume
    c=c+1
fig = plt.figure()
t = np.arange(t_start,t_end,t_step)
plt.plot(t,rho_mean[:c], label='mean')
plt.plot(t,rho_left[:c],label='left')
plt.plot(t,rho_middle[:c],label='middle')
plt.plot(t,rho_right[:c],label='right')
plt.legend()
plt.savefig('mapping_density_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
