#!/usr/bin/env python3

"""
Used to create Couette flow profile plots (average X-vel over simulation time) for a time-parallel simulation.
Plots several PinT iterations in a single plot, can be used to show convergence.

Usage: MaMiCo/build> ../examples/plotting/plot-pint-LBCouette.py

Assumptions:
- this file is in examples/plotting
- there is VTK output from a PinT Couette flow scenario in build
- All *.vtk files in build match with the naming pattern
	LBCouette_supervising_r<RANK>_c<CYCLE>_i<ITERATION>.vtk for G or
	LBCouette_r<RANK>_c<CYCLE>_i<ITERATION>.vtk for F
- No assumptions about total number of ranks, cycles or iterations
- No assumptions about Couette config (wall velocity, viscosity, etc.)

Output will be written to current working directory. Output file name contains total number of solver modes, total number
of MPI ranks, total number of coupling cycles, total number of PinT iterations in data plotted. 
"""
__author__ = "Piet"

import os
import glob
import matplotlib.pyplot as plt
import meshio
import numpy as np
import pandas as pd

def get_attrs(file):
	''' Tries to obtain solver mode, MPI rank, coupling cycle and PinT Iteration from filename '''
	filename = os.path.basename(file)
	parts = filename[:-4].split('_')
	idx = 1
	supervising = (parts[idx] == 'supervising')
	if(supervising):
		idx += 1	
	if parts[idx][0] != 'r':
		raise IOError("Filename error, found strange file " + filename)
	rank = int(parts[idx][1:])
	idx += 1
	if parts[idx][0] != 'c':
		raise IOError("Filename error, found strange file " + filename)
	cycle = int(parts[idx][1:])
	idx += 1
	if(len(parts) <= idx):
		iteration = -1
	else:
		if parts[idx][0] != 'i':
			raise IOError("Filename error, found strange file " + filename)
		iteration = int(parts[idx][1:])
	return (supervising, rank, cycle, iteration)

def get_avg_vel(file):
	''' returns the mean X velocity component from a VTK file'''
	mesh = meshio.read(file)
	avg_vel = np.mean(mesh.cell_data['velocity'], axis=1)
	magnitude = np.linalg.norm(avg_vel)
	return magnitude

def plot(df, label, ls, lw):
	''' plots a single subdataframe '''
	plt.plot(df.cycle, df.avg_vel, label=label, linestyle=ls, linewidth=lw)

def glob_files():
	''' searches for all *.vtk files in ../../build'''
	my_path = os.path.dirname(__file__)
	buildpath = os.path.join(my_path,"../../build")
	filestem = "LBCouette"

	globexpr = os.path.join(buildpath,filestem)  + "*.vtk"
	vtkFilesList = glob.glob(globexpr)
	numfiles = len(vtkFilesList)
	print("Found " + str(numfiles) + " vtk output files in build directory.")
	return vtkFilesList

def get_metadata(files):
	modes = set()
	ranks = set()
	cycles = set()
	iterations = set()
	for file in files:
		(supervising, rank, cycle, iteration) = get_attrs(file)
		modes.add(supervising)
		ranks.add(rank)
		cycles.add(cycle)
		iterations.add(iteration)
	print("Found " + str(len(modes)) + " solver modes, " + str(len(ranks)) + " ranks, "
		+ str(len(cycles)) + " cycles, " + str(len(iterations)) + " iterations.")
	return {'modes': modes, 'ranks': ranks, 'cycles':cycles, 'iterations':iterations}

def read_data(files):
	data = pd.DataFrame(columns=['supervising', 'rank', 'cycle', 'iteration', 'avg_vel'])
	idx = 0
	for file in files:
		(supervising, rank, cycle, iteration) = get_attrs(file)
		vel = get_avg_vel(file)
		data.loc[idx] = [supervising, rank, cycle, iteration, vel]
		idx += 1
	return data.sort_values(by=["supervising", "cycle", "iteration"])

def plot_all(metadata, data):
	for mode in metadata['modes']:
		if mode:
			label = 'G'
			ls = 'dotted'
		else:
			label = 'F'
			ls = 'dashed'
		for it in metadata['iterations']:
			if it >= 0:
				full_label = label + '_' + str(it)
				lw = 1
			else:
				full_label = label + '_sequential'
				lw = 2
				ls = "solid"
			subdataframe = data[(data['supervising'] == mode) & (data['iteration'] == it)]
			size = len(subdataframe.index)
			if(size > 0):
				print("Data series " + full_label + " with size " + str(size) + " is: ")
				print(subdataframe)
				print()
				plot(subdataframe, full_label, ls, lw)
	plt.legend()
	plt.xlabel("Coupling cycle")
	plt.ylabel("Average velocity")
	filename = ('pint' + '_M' + str(len(metadata['modes'])) + '_R' + str(len(metadata['ranks'])) + 
		'_C' + str(len(metadata['cycles'])) + '_I' + str(len(metadata['iterations'])) + '.png')
	plt.savefig(filename)
	plt.show()

def main():
	files = glob_files()
	metadata = get_metadata(files)
	data = read_data(files)
	plot_all(metadata, data)

#only if this file is executed directly i.e.not imported as a module,
#then actually run the simulation
if __name__ == '__main__':
     main()
