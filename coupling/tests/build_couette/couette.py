#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def couette_analytic(x_c, t_c):
    u_w = 0.5
    H = 50.0
    visc = 2.641975309
    sum = 0
    for step in np.arange(1,100):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

x_step = 1.25
colors = [ cm.jet(k) for k in  np.arange(0, 10, 1)]
cmap = plt.get_cmap('jet_r')
x = np.arange(-x_step/2,53,x_step)
for number in np.arange(0, 10, 1):
    f = np.loadtxt("velocity_"+str(4*number)+".txt")
    #print(f)
    plt.plot(x, f, linestyle ='none', marker='x', color=cmap(float(number)/10), label=str(number))
    plt.plot(x, couette_analytic(x,number), color=cmap(float(number)/10))

plt.legend()
plt.savefig('results.png')
plt.show()
