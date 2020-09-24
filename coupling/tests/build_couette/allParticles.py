#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

number = 98

x = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(0), dtype="int")*5+20-2.5
y = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(1), dtype="int")*5+20-2.5
z = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(2), dtype="int")*5+5-2.5
u = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(3), dtype="float")
v = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(4), dtype="float")
w = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(5), dtype="float")
rho = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(6), dtype="float")
fig = plt.figure()
ax = plt.axes(projection="3d")
p=ax.scatter3D(x,y,z,c=w,cmap='Spectral', alpha=0.8)
#p=ax.scatter3D(x,y,z,c=rho/125,cmap='Spectral', alpha=0.8)
print(np.mean(rho/125))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
fig.colorbar(p)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
