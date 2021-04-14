#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

t_start = 0
t_end = 1000
t_step = 2
volume = 5*5*5
ncpd = 12 # number cells per direction
ncpd2 = ncpd*ncpd
ncpd_3 = ncpd-3
iV = 1/volume
datapoints = np.int(t_end/t_step)
rho_mean = np.zeros(datapoints)
rho_left = np.zeros(datapoints)
#rho_left2 = np.zeros(datapoints)
rho_buttom = np.zeros(datapoints)
rho_buttom2 = np.zeros(datapoints)
rho_front = np.zeros(datapoints)
rho_front2 = np.zeros(datapoints)
rho_middle = np.zeros(datapoints)
rho_back = np.zeros(datapoints)
rho_back2 = np.zeros(datapoints)
rho_top = np.zeros(datapoints)
rho_top2 = np.zeros(datapoints)
rho_right = np.zeros(datapoints)
#rho_right2 = np.zeros(datapoints)
c = 0
for t in np.arange(t_start, t_end, t_step):
    rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(t)+".csv", delimiter=';', usecols=(6), dtype="float")
    rho_mean[c] = np.mean(rho)*iV
    rho_left[c] = np.mean(rho[np.arange(ncpd2*3+ncpd*3,ncpd2*3+ncpd*ncpd_3,ncpd)])*iV # 1,4:npcd-3,4
    #rho_left2[c] = np.mean(rho[np.arange(ncpd2*3+ncpd*3+1,ncpd2*3+ncpd*ncpd_3+1,ncpd)])*iV # 2, 4:npcd-3,4
    rho_buttom[c] = np.mean(rho[np.arange(ncpd*3+3,ncpd*ncpd_3+3,ncpd)])*iV # 4,4:npcd-3,1
    rho_buttom2[c] = np.mean(rho[np.arange(ncpd2+ncpd*3+3,ncpd2+ncpd*ncpd_3+3,ncpd)])*iV # 4,4:npcd-3,2
    rho_front[c] = np.mean(rho[np.arange(ncpd2*3+3,ncpd2*3+ncpd_3)])*iV # 4:npcd-3,1,4
    rho_front2[c] = np.mean(rho[np.arange(ncpd2*3+ncpd+3,ncpd2*3+npcd+ncpd_3)])*iV # 4:npcd-3,2,4
    rho_middle[c] = np.mean(rho[np.arange(ncpd2*5+ncpd*3+5,ncpd2*5+ncpd*ncpd-3+5,ncpd)])*iV # 6,4:ncpd-3,6
    rho_back[c] = np.mean(rho[np.arange(ncpd2*3+ncpd2-ncpd+3,ncpd2*3+ncpd2-ncpd+ncpd_3)])*iV # 4:npcd-3,ncpd,4
    rho_back2[c] = np.mean(rho[np.arange(ncpd2*3+ncpd*(ncpd-2)+3,ncpd2*3+ncpd*(ncpd-2)+ncpd_3)])*iV # 4:npcd-3,ncpd-1,4
    rho_top[c] = np.mean(rho[np.arange(ncpd2*(ncpd-1)+ncpd*3+3,ncpd2*(ncpd-1)+ncpd*ncpd_3+3,ncpd)])*iV # 4,4:ncpd-3,12
    rho_top2[c] = np.mean(rho[np.arange(ncpd2*(ncpd-2)+ncpd*3+3,ncpd2*(ncpd-2)+ncpd*ncpd_3+3,ncpd)])*iV #  4,4:ncpd-3,11
    rho_right[c] = np.mean(rho[np.arange(ncpd2*4+ncpd*3+ncpd-1,ncpd2*4+ncpd*ncpd_3+ncpd-1,ncpd)])*iV # ncpd,4:ncpd-3,4
    #rho_right2[c] = np.mean(rho[np.arange(ncpd2*4+ncpd*3+ncpd-2,ncpd2*4+ncpd*ncpd_3+ncpd-2,ncpd)])*iV # ncpd-1,4:ncpd-3,4
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
