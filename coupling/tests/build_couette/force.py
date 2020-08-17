#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

x = np.genfromtxt("force.txt", delimiter=',', usecols=(0), dtype="float")
y = np.genfromtxt("force.txt", delimiter=',', usecols=(1), dtype="float")
z = np.genfromtxt("force.txt", delimiter=',', usecols=(2), dtype="float")
f = np.genfromtxt("force.txt", delimiter=',', usecols=(3), dtype="float")
fig = plt.figure()
ax = plt.axes(projection="3d")
p=ax.scatter3D(x[500:1000],y[500:1000],z[500:1000],c=f[500:1000],cmap='Greens', alpha=0.8)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
fig.colorbar(p)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
