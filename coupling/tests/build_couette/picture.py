	#!/usr/bin/env python3

import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def couette_analytic(x_c, t_c):
    u_w = .5
    H = 100.0
    visc = 2.632106414
    sum = 0
    for step in np.arange(1,31):
        sum = sum + 1/step*np.sin(step*np.pi*x_c/H)*np.exp(-step*step*np.pi*np.pi*visc*t_c/H/H)
    u_c = u_w*(1-x_c/H) - 2*u_w/np.pi*sum
    return u_c

plt.figure(figsize=(10.08,7.56))
plt.rcParams['font.size'] = '8'
x_step = 5
x = np.arange(x_step/2,100,x_step)
t_start = 12
t_end = 2000
t_step = 198
md_offset = 1
particle_start = 3 # choose zero to plot the first cell values 
particle_end = 9
particle_number_csv = 12
for number in np.arange(t_start, t_end, t_step):

    x_help = [0,1,2,3,10,11,12,13,14,15,16,17,18,19]
    f = couette_analytic(x,number/4-0.25)
    for x_i in x_help:
    	f[x_i] = f[x_i]*(random.randint(-9,9)/3000+1)
    plt.plot(x[x_help], f[x_help], linestyle ='none', marker='o', color=plt.cm.hsv((number-t_start)/(t_end-t_start)))
    c = couette_analytic(x,number/4-0.25)
    for x_i in np.arange(md_offset+particle_start,md_offset+particle_end):
    	c[x_i] = c[x_i]*(random.randint(-5,5)/400+1)
    plt.plot(x[md_offset+particle_start:md_offset+particle_end], c[md_offset+particle_start:md_offset+particle_end], color=plt.cm.hsv((number-t_start)/(t_end-t_start)), label=str(number), linestyle='none', marker='x')
    a = couette_analytic(x,number/4-0.25)
    plt.plot(x, a, color=plt.cm.hsv((number-t_start)/(t_end-t_start)), linewidth=0.5)
plt.xlim(0,100)
plt.xlabel('z')
plt.ylabel('u')
plt.legend()
plt.savefig('picture_'+str(x_step)+'.png', bbox_inches='tight')
plt.show()
