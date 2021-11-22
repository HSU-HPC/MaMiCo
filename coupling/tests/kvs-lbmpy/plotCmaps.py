#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs 
from matplotlib import rcParams
from pandas import read_csv

fs = 14
rcParams.update({'font.size': fs, 'font.weight': 'bold'})

fig, ax = plt.subplots(1, 3)
fig.set_size_inches(12,5)

# 'noisy_signal', 
files = ['original_signal', 'NLM_old', 'NLM_new']
data = [read_csv(file+".csv", delimiter=";", usecols=[0,1,2,3,4], \
	names=["t","x","y","z","m"], index_col=None) for file in files]

# 'Noisy Data', 
titles = ['Original Signal', 'NLM (1)', 'NLM (2)']

#Plot all pcolor plots
for i in range(3):
	timeslice = data[i].loc[data[i]['t'].isin(range(50,58))]
	xslice = timeslice.loc[(timeslice['y'] == 5) & (timeslice['z'] == 5)]
	a = np.array(xslice.loc[xslice['x'].isin(range(4,10))]['m'])
	a.shape = (a.size//6, 6)
	pc = ax[i].pcolor(a, vmin=0.85, vmax=1.15)
	ax[i].set_title(titles[i],fontsize=fs, fontweight='bold')
	ax[i].set_xlabel('cell x',fontsize=fs, fontweight='bold')
	ax[i].set_ylabel('time t (coupling cycles)',fontsize=fs, fontweight='bold')
	ax[i].set_xticks([0.5,1.5,2.5,3.5, 4.5, 5.5])
	ax[i].set_xticklabels(['0','1','2','3','4','5'],fontsize=fs, fontweight='bold')
	ax[i].set_yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5])
	ax[i].set_yticklabels(['0','1','2','3','4','5','6','7'],fontsize=fs, fontweight='bold')

fig.subplots_adjust(right=0.88, wspace=0.3)
#cbar_ax = fig.add_axes([0.90, 0.15, 0.01, 0.7]) #left, bottom, width, height
#cbar = fig.colorbar(pc, cax=cbar_ax)
#cbar.set_label("density",fontsize=fs, fontweight='bold')

#fig.tight_layout()
plt.show()
