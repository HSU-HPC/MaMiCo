#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t_start = 0
t_end = 1000
t_step = 2
volume = 5*5*5
rho_mean = np.zeros(np.int(t_end/t_step))
rho_left = np.zeros(np.int(t_end/t_step))
#rho_left2 = np.zeros(np.int(t_end/t_step))
rho_buttom = np.zeros(np.int(t_end/t_step))
rho_buttom2 = np.zeros(np.int(t_end/t_step))
rho_front = np.zeros(np.int(t_end/t_step))
rho_front2 = np.zeros(np.int(t_end/t_step))
rho_middle = np.zeros(np.int(t_end/t_step))
rho_back = np.zeros(np.int(t_end/t_step))
rho_back2 = np.zeros(np.int(t_end/t_step))
rho_top = np.zeros(np.int(t_end/t_step))
rho_top2 = np.zeros(np.int(t_end/t_step))
rho_right = np.zeros(np.int(t_end/t_step))
#rho_right2 = np.zeros(np.int(t_end/t_step))
c = 0
for t in np.arange(t_start, t_end, t_step):
    rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(t)+".csv", delimiter=';', usecols=(6), dtype="float")
    rho_mean[c] = np.mean(rho/volume)
    rho_left[c] = np.mean(rho[np.arange(744,829,12)])/volume
    #rho_left2[c] = np.mean(rho[np.arange(745,830,12)])/volume
    rho_buttom[c] = np.mean(rho[np.arange(29,114,12)])/volume
    rho_buttom2[c] = np.mean(rho[np.arange(173,258,12)])/volume
    rho_front[c] = np.mean(rho[np.arange(581,1014,144)])/volume
    rho_front2[c] = np.mean(rho[np.arange(593,1026,144)])/volume
    rho_middle[c] = np.mean(rho[np.arange(749,834,12)])/volume
    rho_back[c] = np.mean(rho[np.arange(713,1146,144)])/volume
    rho_back2[c] = np.mean(rho[np.arange(701,1134,144)])/volume
    rho_top[c] = np.mean(rho[np.arange(1613,1698,12)])/volume
    rho_top2[c] = np.mean(rho[np.arange(1469,1554,12)])/volume
    rho_right[c] = np.mean(rho[np.arange(755,840,12)])/volume
    #rho_right2[c] = np.mean(rho[np.arange(756,841,12)])/volume
    c=c+1
fig = plt.figure()
t = np.arange(t_start,t_end,t_step)
plt.plot(t,rho_mean[:c], label='mean')
#plt.plot(t,rho_left[:c],label='left')
#plt.plot(t,rho_left2[:c],label='left 2')
plt.plot(t,rho_buttom[:c],label='buttom')
plt.plot(t,rho_buttom2[:c],label='buttom 2')
#plt.plot(t,rho_front[:c],label='front')
#plt.plot(t,rho_front2[:c],label='front 2')
plt.plot(t,rho_middle[:c],label='middle')
#plt.plot(t,rho_back[:c],label='back')
#plt.plot(t,rho_back2[:c],label='back 2')
plt.plot(t,rho_top[:c],label='top')
plt.plot(t,rho_top2[:c],label='top 2')
#plt.plot(t,rho_right[:c],label='right')
#plt.plot(t,rho_right2[:c],label='right 2')
plt.legend()
plt.savefig('mapping_density_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
