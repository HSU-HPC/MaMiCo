#!/usr/bin/env python3

import numpy as np
from pandas import read_csv
import math

original_file = 'lbm'
filtered_files = [ "Gaussian-0.5", "Gaussian-1", "Gaussian-1.5", "Gaussian-2", "Gaussian-2.5", "Gaussian-3", "Gaussian-3.5", "Gaussian-4", "Gaussian-4.5", "Gaussian-5", \
	'POD-5-1', 'POD-5-2', 'POD-5-3', 'POD-20-1', 'POD-20-2', 'POD-20-3', 'POD-35-1', 'POD-35-2', 'POD-35-3', 'POD-50-1', 'POD-50-2', 'POD-50-3', 'POD-65-1', 'POD-65-2', 'POD-65-3', 'POD-80-1', 'POD-80-2', 'POD-80-3', \
	'NLM_0.1', 'NLM_1', 'NLM_10', 'NLM_100', 'NLM_1000', 'NLM_10000']

#TODO: what even are the columns in lbm.csv?
original_data = read_csv(original_file+".csv", delimiter=";", usecols=[0,1], names=["t", "m"], index_col=None)
filtered_data = [read_csv(file+".csv", delimiter=";", usecols=[0,1,2,3,7], \
	names=["t","x","y","z","m"], index_col=None) for file in filtered_files]


for i in range(len(filtered_data)):
	rows = filtered_data[i]#.loc[data[i]['t'].isin(range(50,100))]

	#select cell 129
	rows = rows.loc[rows['x'] == 3]
	rows = rows.loc[rows['y'] == 3]
	rows = rows.loc[rows['z'] == 3]

	filtered_masses = np.array(rows['m'])

	diff = np.array(original_data['m']) - filtered_masses
	MSE = np.dot(diff, diff) / filtered_masses.size

	print(filtered_files[i] + ": MSE = " + str(MSE))
