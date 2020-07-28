#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def couette_analytic(x_c, t_c):
    u_w = 0.5
    H = 50.0
    visc = 2.641975309
    sum = 0
    for step in np.arange(1,31):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

x_step = 2.5
x = np.arange(x_step/2,50,x_step)
t_start = 1
t_end = 11
t_step = 1
particle_start = 3
particle_number = 6
for number in np.arange(t_start, t_end, t_step):
    c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(3), dtype="float")
    c_0 = np.zeros(particle_number)
    help = particle_number*particle_number
    for i in np.arange(0, particle_number):
        c_0[i] = np.mean(c[np.int(i*help):np.int((i+1)*help)])
    f = np.loadtxt("velocity_"+str(number)+".txt")
    plt.plot(x, f, linestyle ='none', marker='x', color=plt.cm.hsv((number-t_start)/(t_end-t_start)), label=str(number))
    #c = couette_analytic(x,number/4)
    plt.plot(x[particle_start:particle_start+particle_number], c_0, color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
    #plt.plot(x, c, color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
    #print(couette_analytic(x,10/4))
plt.legend()
plt.savefig('particle_'+str(x_step)+'.png')
plt.show()
