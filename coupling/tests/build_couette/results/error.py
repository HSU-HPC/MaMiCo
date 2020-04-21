#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

g = np.loadtxt("error_25.txt")
plt.plot(g[10:], color='black', label='2.5')
f = np.loadtxt("error_4.txt")
plt.plot(f[10:], color='blue', label='4')
j = np.loadtxt("error_2.txt")
plt.plot(j[10:], color='red', label='2')
h = np.loadtxt("error_5.txt")
plt.plot(h[10:], color='green', label='5')
#print(h/g)
plt.ylabel('max error')
plt.xlabel('timestep')
plt.legend()
plt.savefig('max_errror.png', bbox_inches='tight')
plt.show()
