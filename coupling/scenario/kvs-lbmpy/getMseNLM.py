#!/usr/bin/env python3

import math

import numpy as np
from pandas import read_csv

files = ['original_signal.105', 'NLM_junction']
data = [read_csv(file+".csv", delimiter=";", usecols=[0,1,2,3,7], \
	names=["t","x","y","z","m"], index_col=None) for file in files]

a = [0,0]
for i in range(2):
	rows = data[i].loc[data[i]['t'].isin(range(50,100))]

	rows = rows.loc[rows['x'].isin(range(1,5))]
	rows = rows.loc[rows['y'].isin(range(1,5))]
	rows = rows.loc[rows['z'].isin(range(1,5))]
	a[i] = np.array(rows['m'])

diff = a[0]-a[1]
MSE = np.dot(diff, diff) / a[0].size

print("MSE = " + str(MSE))
