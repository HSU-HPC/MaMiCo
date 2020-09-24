#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def couette_analytic(x_c, t_c):
    u_w = 1.5
    H = 100.0
    visc = 2.632106414
    sum = 0
    for step in np.arange(1,31):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

plt.figure(figsize=(12,9))
x_step = 5
x = np.arange(x_step/2,100,x_step)
t_start = 4
t_end = 98
t_step = 10
particle_start = 1
particle_number = 12
for number in np.arange(t_start, t_end, t_step):
    f = np.genfromtxt("velocity_"+str(number)+".txt", usecols=(0), delimiter=',', dtype="float")
    x_help = [0,1,2,3,10,11,12,13,14,15,16,17,18,19]
    plt.plot(x[x_help], f, linestyle ='none', marker='o', color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
    x_c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(0), dtype="int")
    y_c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(1), dtype="int")
    z_c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(2), dtype="int")
    c = np.genfromtxt("CouetteAvgMultiMDCells_0_"+str(number)+".csv", delimiter=';', usecols=(3), dtype="float")
    c_0 = np.zeros(particle_number)
    c_c = np.zeros(particle_number)
    help = particle_number*particle_number
    for i in np.arange(0, particle_number):
        for k in np.arange(i*help,(i+1)*help):
            #if np.logical_and(not x_c[k]==1,not x_c[k]==12):
                #if np.logical_and(not y_c[k]==1,not y_c[k]==12):
            c_0[i] = c_0[i]+c[k]
            c_c[i] = c_c[i]+1
    c_0=c_0/c_c
    plt.plot(x[particle_start:particle_start+particle_number], c_0, color=plt.cm.hsv((number-t_start)/(t_end-t_start)), label=str(number), linestyle=':', marker='x')
    # help = np.arange(65, 1652, 144)
    # plt.plot(x[np.arange(1,13)], c[help], color=plt.cm.hsv((number-t_start)/(t_end-t_start)), label=str(number), linestyle=':', marker='x')
    a = couette_analytic(x,number/4-0.25)
    plt.plot(x, a, color=plt.cm.hsv((number-t_start)/(t_end-t_start)), linewidth=0.5)
plt.legend()
plt.savefig('particle_'+str(x_step)+'.png', bbox_inches='tight')
plt.show()
