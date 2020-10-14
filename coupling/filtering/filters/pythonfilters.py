#!/usr/bin/env python3

# This file is part of the Mamico project
# Author: Felix Maurer
# September 2020
# Helmut Schmidt University, Hamburg. Chair for High Performance Computing
# BSD license, see the copyright notice in Mamico's main folder

import numpy as np
from scipy import fftpack
#TODO: Clashes with other plots
import matplotlib.pyplot as plt

# Collection of filter classes written in Python. 
# These cannot be added via .xml configuration but rather via mamico.MacrosocopicCellService.addFilterToSequence() at runtime.
# Can only be used from within a Python context.

DEBUG_STROUHAL_PYTHON = True

#Global filter functions that can find application regardless of filter in use:
def returnCellData(cellData):
    print("Copying cell data.")
    return cellData


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
    def addDataPoint(self, cellData, indices):
        total = 0
        for cell in cellData:
            total = total + cell[1]
        self.avgs.append(total/len(cellData))
     
        if DEBUG_STROUHAL_PYTHON:
            print("STROUHAL-PYTHON: New datapoint added: " + str(total/len(cellData)))

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
