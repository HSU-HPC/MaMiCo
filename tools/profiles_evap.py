#!/usr/bin/env python3
import sys, getopt
import numpy as np
import pandas as pd
import matplotlib as mpl
import os
from matplotlib import pyplot as plt
import seaborn as sns

# simulation parameters
dt = 0.1
max_cc = 150
ccs = np.arange(1,max_cc)

title = 'Velocity profiles'
filename = 'profiles.eps'

def load_mamico_data(csv_file):
	""" Get CSV data from one cycle for mamico
	"""
	# load data in pandas DataFrame
	df = pd.read_csv(csv_file, sep=";")
	return df

def plot_one_coupling_cycle(cc, color, ax, path):
	""" Plot MaMiCo mass
	"""
	csv_file = os.path.join(path, 'mamico', 'results', 'results_0_'+str(cc)+'.csv')
	res = load_mamico_data(csv_file)
	ax.plot(cc, res.loc[:,'m_micro'].sum(), color=color, marker='.', linestyle='none', markersize=4)
	""" Mass
	"""
	#print('z:'+str(res.loc[:,'z']))
	#print('rho:'+str(res.loc[:,'avg_mass'].sum()/volume_inner_md_domain))

def main(argv):
	opts, args = getopt.getopt(argv,'hp::',['path='])
	path = '.'
	for opt, arg in opts:
		if opt == '-h':
			print('profiles.py -p path')
			sys.exit()
		elif opt in ('-p', '--path'):
			path = arg
	fig = plt.figure(figsize=(7,5))
	ax = fig.add_subplot(111)
	cmap = mpl.colormaps['Spectral']
	norm = mpl.colors.Normalize(vmin=ccs[0], vmax=ccs[-1])

	for cc in ccs:
		plot_one_coupling_cycle(cc, cmap(norm(cc)), ax, path)

	ax.set_xlim(0,max_cc)
	ax.set_xlabel('cc', fontsize=14)
	ax.set_ylabel('microscopic mass', fontsize=14)
	# ax.set_title(title, fontsize=14)
	ax.grid(visible=True)
	ax.legend()
	plt.show()
	# fig.savefig(filename, dpi=250)

if __name__ == "__main__":
   main(sys.argv[1:])
