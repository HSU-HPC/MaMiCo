#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

number = 570

numberCells = 14
divFactor = np.ones([3]);
divFactor[1] = divFactor[0]*numberCells
divFactor[2] = divFactor[1]*numberCells

def index2vector(index2):
    vector = np.zeros([np.size(index2),3])
    for k in np.arange(np.size(index2)):
        help = index2[k]
        for i in [2,1]:
            vector[k,i] = np.int(help/divFactor[i])
            help = np.int(help - vector[k,i]*divFactor[i])
        vector[k,0] = help
    return vector

index = np.genfromtxt("excessMass_"+str(number)+".txt", delimiter=';', usecols=(0), dtype="int")
mass = np.zeros([np.size(index),3])
for i in np.arange(3):
    mass[:,i] = np.genfromtxt("excessMass_"+str(number)+".txt", delimiter=';', usecols=(1), dtype="float")
vectorIndex=index2vector(index)
fig = plt.figure()
ax = plt.axes(projection="3d")
p=ax.scatter3D(vectorIndex[:,0],vectorIndex[:,1],vectorIndex[:,2],c=np.mean(mass,axis=1),cmap='RdBu', alpha=0.8,vmin=-3, vmax=3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
fig.colorbar(p)
#plt.savefig('allParticle_'+str(t_end)+'.png', bbox_inches='tight')
plt.show()
