#!/usr/bin/env python3

import os
import glob

my_path = os.path.dirname(__file__)
buildpath = os.path.join(my_path,"../../build")
filestem = "LBCouette"

globexpr = os.path.join(buildpath,filestem)  + "*.vtk"
vtkFilesList = glob.glob(globexpr)
numfiles = len(vtkFilesList)

print("Found " + str(numfiles) + " vtk output files in build directory.")

def get_attrs(file):
	filename = os.path.basename(file)
	parts = filename[:-4].split('_')
	idx = 1
	supervising = (parts[idx] == 'supervising')
	if(supervising):
		idx += 1	
	if parts[idx][0] != 'r':
		raise IOException("Filename error, found strange file " + filename)
	rank = int(parts[idx][1:])
	idx += 1
	if parts[idx][0] != 'c':
		raise IOException("Filename error, found strange file " + filename)
	cycle = int(parts[idx][1:])
	idx += 1
	if(len(parts) <= idx):
		iteration = 0
	else:
		if parts[idx][0] != 'i':
			raise IOException("Filename error, found strange file " + filename)
	iteration = int(parts[idx][1:])
	return (supervising, rank, cycle, iteration)

modes = set()
ranks = set()
cycles = set()
iterations = set()
for file in vtkFilesList:
	(supervising, rank, cycle, iteration) = get_attrs(file)
	modes.add(supervising)
	ranks.add(rank)
	cycles.add(cycle)
	iterations.add(iteration)

print("Found " + str(len(modes)) + " solver modes, " + str(len(ranks)) + " ranks, "
	+ str(len(cycles)) + " cycles, " + str(len(iterations)) + " iterations.")

import meshio
import numpy as np
def get_avg_vel(file):
	mesh = meshio.read(file)
	avg_vel = np.mean(mesh.cell_data['velocity'], axis=1)
	magnitude = np.linalg.norm(avg_vel)
	return magnitude

#for file in vtkFilesList:
#	print(os.path.basename(file) + " = " + str(get_avg_vel(file)))

import pandas as pd
data = pd.DataFrame(columns=['supervising', 'rank', 'cycle', 'iteration', 'avg_vel'])
idx = 0
for file in vtkFilesList:
	(supervising, rank, cycle, iteration) = get_attrs(file)
	vel = get_avg_vel(file)
	data.loc[idx] = [supervising, rank, cycle, iteration, vel]
	idx += 1
data = data.sort_values(by=["supervising", "cycle", "iteration"])
#print(data)

import matplotlib.pyplot as plt
import random

def plot(df, label):
	lw=random.randint(1, 5)
	ls=random.choice(['-','--','-.',':'])
	plt.plot(df.cycle, df.avg_vel, linestyle=ls, linewidth=lw, label=label)

for mode in modes:
	if mode:
		label = 'G'
	else:
		label = 'F'
	for it in iterations:
		full_label = label + '_' + str(it)
		subdataframe = data[(data['supervising'] == mode) & (data['iteration'] == it)]
		size = len(subdataframe.index)
		if(size > 0):
			print("Data series " + full_label + " with size " + str(size) + " is: ")
			print(subdataframe)
			print()
			plot(subdataframe, full_label)

plt.legend()
plt.xlabel("Coupling cycle")
plt.ylabel("Average velocity")
plt.savefig('pint.png')
plt.show()
