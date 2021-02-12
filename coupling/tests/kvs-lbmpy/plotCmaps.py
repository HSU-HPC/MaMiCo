#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs 
from pandas import read_csv

fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(1, 5)

dataCS = np.random.rand(8, 4) #TODO
dataMD = read_csv("testPlotCmaps.csv", delimiter=";", usecols=[0,1,7], names=["Iteration", "X", "Mass"], index_col=None)
dataPOD = read_csv("testPlotCmaps.csv", delimiter=";", usecols=[0,1,7], names=["Iteration", "X", "Mass"], index_col=None)
dataGAUSS = read_csv("testPlotCmaps.csv", delimiter=";", usecols=[0,1,7], names=["Iteration", "X", "Mass"], index_col=None)
dataNLM = read_csv("testPlotCmaps.csv", delimiter=";", usecols=[0,1,7], names=["Iteration", "X", "Mass"], index_col=None)

#Plot all pcolor plots

#CS
ax0.pcolor(dataCS)
ax0.set_title('Macro Solver')

#MD
a = np.array(dataMD["Mass"])
a.shape = (a.size//4, 4)
ax1.pcolor(a)
ax1.set_title('MD')

#POD
a = np.array(dataPOD["Mass"])
a.shape = (a.size//4, 4)
ax2.pcolor(a)
ax2.set_title('POD')

#GAUSS
a = np.array(dataGAUSS["Mass"])
a.shape = (a.size//4, 4)
ax3.pcolor(a)
ax3.set_title('Gaussian')

#NLM
a = np.array(dataNLM["Mass"])
a.shape = (a.size//4, 4)
ax4.pcolor(a)
ax4.set_title('NLM')



fig.tight_layout()
plt.show()
