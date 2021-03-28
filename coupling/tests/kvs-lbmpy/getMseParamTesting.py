#!/usr/bin/env python3

import numpy as np
from pandas import read_csv
import math
import matplotlib.pyplot as mplt

original_file = 'lbm'
filtered_files = [ ["Gaussian-0.5", "Gaussian-1", "Gaussian-1.5", "Gaussian-2", "Gaussian-2.5", "Gaussian-3", "Gaussian-3.5", "Gaussian-4", "Gaussian-4.5", "Gaussian-5"], \
	['POD-20-1', 'POD-35-1', 'POD-50-1', 'POD-65-1', 'POD-80-1'],\
	['POD-20-2', 'POD-35-2', 'POD-50-2', 'POD-65-2', 'POD-80-2'],\
	['POD-20-3', 'POD-35-3', 'POD-50-3', 'POD-65-3', 'POD-80-3'],\
	['NLM_0.1', 'NLM_1', 'NLM_10', 'NLM_100'] ]

original_data = read_csv(original_file+".csv", delimiter=";", usecols=[0,1,2], names=["t", "x", "y"], index_col=None)

MSE_x = []
MSE_y = []

#for each filter type
for f in filtered_files:
	filter_data = [read_csv(file+".csv", delimiter=";", usecols=[0,1,2,3,7,8,9,10], \
	names=["t","x","y","z","m", "px", "py", "pz"], index_col=None) for file in f]

	#new MSE lists for this filter type
	MSE_x.append([])
	MSE_y.append([])
	print(MSE_x)
		
	#for each instance of that filter type
	for i in range(len(filter_data)):


		rows = filter_data[i]#.loc[data[i]['t'].isin(range(50,100))]

		#select cell 129
		rows = rows.loc[rows['x'] == 3]
		rows = rows.loc[rows['y'] == 3]
		rows = rows.loc[rows['z'] == 3]

		filtered_px = np.array(rows['px'])
		filtered_py = np.array(rows['py'])
		filtered_m = np.array(rows['m'])

		#compute velocity
		filtered_velx = filtered_px / filtered_m
		filtered_vely = filtered_py / filtered_m

		#compute MSE in x,y dirs
		diff_x = np.array(original_data["x"]) - filtered_velx
		diff_y = np.array(original_data["y"]) - filtered_vely

		MSE_x[-1].append(np.dot(diff_x, diff_x) / diff_x.size)
		MSE_y[-1].append(np.dot(diff_y, diff_y) / diff_y.size)

		print(f[i] + ": MSE (vel x) = " + str(MSE_x[-1][-1]))
		print(f[i] + ": MSE (vel y) = " + str(MSE_y[-1][-1]))
	
def plot_mses(file_location, tested_parameter, x_axis, f_index):
	mplt.style.use("seaborn")
	print(MSE_y[f_index])
	fig, ax = mplt.subplots(2)
	ax[0].plot(x_axis, MSE_x[f_index], "-", label = file_location)
	ax[1].plot(x_axis, MSE_y[f_index], "-", label = file_location)

	#labeling
	ax[0].set_ylabel("MSE of x velocity")
	ax[1].set_ylabel("MSE of y velocity")
	ax[0].set_xlabel(tested_parameter)
	ax[1].set_xlabel(tested_parameter)

	fig.tight_layout()

	file_location = file_location + ".png"
	mplt.savefig(file_location)

plot_mses("Gaussian", "sigma", [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5], 0)
plot_mses("POD-1", "tws", [20, 35, 50, 65, 80], 1)
plot_mses("POD-2", "tws", [20, 35, 50, 65, 80], 2)
plot_mses("POD-3", "tws", [20, 35, 50, 65, 80], 3)
plot_mses("NLM", "sigma²", ["0.1", "1", "10", "100"], 4)
