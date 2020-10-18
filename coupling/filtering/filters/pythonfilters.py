#!/usr/bin/env python3

# This file is part of the Mamico project
# Author: Felix Maurer
# September 2020
# Helmut Schmidt University, Hamburg. Chair for High Performance Computing
# BSD license, see the copyright notice in Mamico's main folder

import numpy as np
#TODO: Clashes with other plots
import matplotlib.pyplot as plt

# Collection of filter classes written in Python. 
# These cannot be added via .xml configuration but rather via mamico.MacrosocopicCellService.addFilterToSequence() at runtime.
# Can only be used from within a Python context.

DEBUG_PYTHON_FILTERS = False

#Global filter functions that can find application regardless of filter in use:
def returnCellData(cellData):
    if DEBUG_PYTHON_FILTERS:
        print("Copying cell data.")

    return cellData


from scipy import fftpack
DEBUG_STROUHAL_PYTHON = True

# WARNING: ONLY USE THIS IN A VELOCITY-ONLY FILTER SEQUENCE.
class StrouhalPython():
    def __init__ (self, u, d):
        self.u = u
        self.d = d
        self.avgs = []
        if DEBUG_STROUHAL_PYTHON:
            print("STROUHAL-PYTHON: New filter created. u/d: " + str(u) + "/" + str(d))

    # If a FilterFromFunction-instance's operator() call calls this more than once (e.g. the sequence filters not exclusively velocity),
    # this method creates faulty data points. => Only use this filter in a velocity-only filter sequence.
    def addDataPoint(self, cellData):
        shape = cellData.shape
        if shape[0] == 0 and shape[1] == 0 and shape[2] == 0:
            print("STROUHAL-PYTHON: WARNING: EMPTY INPUT ARRAY. NO DATAPOINT ADDED.")
            return cellData
        
        #else
        num_cells = shape[0] * shape[1] * shape[2] 
        total = 0
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    #only y value is relevant
                    total = total + cellData[x,y,z,1]
        avg = total/(num_cells)
        self.avgs.append(avg)
        if DEBUG_STROUHAL_PYTHON:
            print("STROUHAL-PYTHON: New datapoint added: " + str(avg))

        return cellData

    def calculateStrouhalNumber(self):
        plt.cla()
        #vels, freqs = plt.subplots()

        x = np.arange(len(self.avgs))
        plt.plot(x, self.avgs)
        
        plt.cla()
        
        f = fftpack.fftfreq(len(self.avgs))
        dft = fftpack.fft(self.avgs)
        plt.plot(f, dft)

        plt.savefig("python_plot.png")
        print("STROUHAL-PYTHON: Done plotting. ")

from scipy.ndimage import gaussian_filter
DEBUG_GAUSS_PYTHON = True

class GaussPython():
    def __init(self, sigma):
        self.sigma = sigma
        if DEBUG_GAUSS_PYTHON:
            print("GAUSS-PYTHON: New filter created. Sigma: " + str(sigma))

    def apply(self, cellData):
        return gaussian_filter(cellData, self.sigma)


