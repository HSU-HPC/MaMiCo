#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def couette_analytic(x_c, t_c):
    u_w = 1.5
    H = 50.0
    visc = 2.641975309
    sum = 0
    for step in np.arange(1,31):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

x_step = 2.5
x = np.arange(x_step/2,50,x_step)
t_start = 15
t_end = 50
t_step = 4
particle_start = 3
particle_number = 6
for number in np.arange(t_start, t_end, t_step):
    c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(3), dtype="float")
    # c_0 = np.zeros(particle_number)
    # help = particle_number*particle_number
    # for i in np.arange(0, particle_number):
    #     c_0[i] = np.mean(c[np.int(i*help):np.int((i+1)*help)])
    help = np.arange(14,216, 36)
    plt.plot(x[particle_start:particle_start+particle_number], c[help], color=plt.cm.hsv((number-t_start)/(t_end-t_start)), label=str(number), linestyle=':', marker='x')
    f = np.genfromtxt("velocity_"+str(number+1)+".txt", usecols=(0), delimiter=',', dtype="float")
    x_help = [0,1,2,9,10,11,12,13,14,15,16,17,18,19]
    plt.plot(x[x_help], f, linestyle ='none', marker='o', color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
    a = couette_analytic(x,number/4)

    plt.plot(x, a, color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
plt.legend()
plt.savefig('particle_'+str(x_step)+'.png')
plt.show()