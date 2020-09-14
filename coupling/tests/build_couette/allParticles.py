#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

number = 40

def couette_analytic(x_c):
    u_w = 1.5
    H = 50.0
    visc = 2.632106414
    t_c = number*0.25
    sum = 0
    for step in np.arange(1,31):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

x = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(0), dtype="int")*2.5+10-1.25
y = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(1), dtype="int")*2.5+10-1.25
z = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(2), dtype="int")*2.5-1.25
u = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(3), dtype="float")
v = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(4), dtype="float")
w = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(5), dtype="float")
rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(6), dtype="float")
fig = plt.figure()
ax = plt.axes(projection="3d")
#p=ax.scatter3D(x,y,z,c=(u-couette_analytic(z))/couette_analytic(z),cmap='Spectral', alpha=0.8, vmin=-5,vmax=5)
#p=ax.scatter3D(x,y,z,c=u,cmap='Spectral', alpha=0.8)
p=ax.scatter3D(x,y,z,c=rho/15.625,cmap='Spectral', alpha=0.8)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
fig.colorbar(p)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
