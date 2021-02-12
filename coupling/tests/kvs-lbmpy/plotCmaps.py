#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs 
from pandas import read_csv

fig, ax = plt.subplots(1, 5)
fig.set_size_inches(12,5)

files = ['original_signal', 'noisy_signal', 'POD', 'Gaussian', 'NLM_junction']
data = [read_csv(file+".csv", delimiter=";", usecols=[0,1,2,3,7], \
	names=["t","x","y","z","m"], index_col=None) for file in files]

titles = ['Original Signal', 'Noisy Data', 'POD', 'Gaussian', 'NLM']

#Plot all pcolor plots
for i in range(5):
	timeslice = data[i].loc[data[i]['t'].isin(range(50,58))]
	xslice = timeslice.loc[(timeslice['y'] == 1) & (timeslice['z'] == 1)]
	a = np.array(xslice.loc[xslice['x'].isin(range(1,5))]['m'])
	a.shape = (a.size//4, 4)
	pc = ax[i].pcolor(a, vmin=0.85, vmax=1.15)
	ax[i].set_title(titles[i])
	ax[i].set_xlabel('cell x')
	ax[i].set_ylabel('time t')

fig.subplots_adjust(right=0.88, wspace=0.32)
cbar_ax = fig.add_axes([0.90, 0.15, 0.01, 0.7]) #left, bottom, width, height
cbar = fig.colorbar(pc, cax=cbar_ax)
cbar.set_label("density")

#fig.tight_layout()
plt.show()
