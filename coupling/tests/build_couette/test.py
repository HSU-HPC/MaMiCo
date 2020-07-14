#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

c = np.loadtxt("test.txt")
divFactor = np.ones([3])
divFactor[1]=14
divFactor[2]=14*14
vector = np.zeros([3,np.size(c)])
for i in np.arange(np.size(c)):
    vector[2,i] = np.trunc(c[i]/divFactor[2])
    help = c[i] - vector[2,i]*divFactor[2]
    vector[1,i] = np.trunc(help/divFactor[1])
    help = help - vector[1,i]*divFactor[1]
    vector[0,i] = help

g = np.loadtxt("tost.txt")
voctor = np.zeros([3,np.size(g)])
for i in np.arange(np.size(g)):
    voctor[2,i] = np.trunc(g[i]/divFactor[2])
    help = g[i] - voctor[2,i]*divFactor[2]
    voctor[1,i] = np.trunc(help/divFactor[1])
    help = help - voctor[1,i]*divFactor[1]
    voctor[0,i] = help


ax = plt.axes(projection='3d')
ax.scatter3D(vector[0], vector[1], vector[2]);
ax.scatter3D(voctor[0], voctor[1], voctor[2], c='red');

plt.show()
