#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

x1 = np.genfromtxt("index_0.txt", delimiter=' ', usecols=(0), dtype="float")
y1 = np.genfromtxt("index_0.txt", delimiter=' ', usecols=(1), dtype="float")
z1 = np.genfromtxt("index_0.txt", delimiter=' ', usecols=(2), dtype="float")

fig = plt.figure()
ax = plt.axes(projection="3d")
#p=ax.plot(x1,y1,z1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

x2 = np.genfromtxt("index_1.txt", delimiter=' ', usecols=(0), dtype="float")
y2 = np.genfromtxt("index_1.txt", delimiter=' ', usecols=(1), dtype="float")
z2 = np.genfromtxt("index_1.txt", delimiter=' ', usecols=(2), dtype="float")
#p=ax.plot(x2,y2,z2)

x3 = np.genfromtxt("index_2.txt", delimiter=' ', usecols=(0), dtype="float")
y3 = np.genfromtxt("index_2.txt", delimiter=' ', usecols=(1), dtype="float")
z3 = np.genfromtxt("index_2.txt", delimiter=' ', usecols=(2), dtype="float")
#p=ax.plot(x3,y3,z3)

x4 = np.genfromtxt("index_3.txt", delimiter=' ', usecols=(0), dtype="float")
y4 = np.genfromtxt("index_3.txt", delimiter=' ', usecols=(1), dtype="float")
z4 = np.genfromtxt("index_3.txt", delimiter=' ', usecols=(2), dtype="float")
p=ax.plot(x4,y4,z4)

plt.show()
