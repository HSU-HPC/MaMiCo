#!/usr/bin/env python3
import sys, getopt
import numpy as np
import pandas as pd
import matplotlib as mpl
import os
from matplotlib import pyplot as plt
import seaborn as sns

## simulation parameters
# macro timestep
dt = 0.25
# coupling cycles to plot
ccs = [100, 200, 500, 1000]
# volume of the inner MD domain
volume_inner_md_domain = 6*6*6*2.5*2.5*2.5
# height
H = 50
# wall velocity
u_w = 0.5
# z cell index of the inner MD domain
index=[4,5,6,7,8,9]
# title of the plot
title = 'Velocity profiles'
# filename for export
filename = 'profiles.eps'
       
def couette_analytic(z, t):
	""" Analytic Couette startup equation
	"""
	# kinematic viscosty = dynamic viscosity / density
	v = 2.14 / .813037037
 
	k_sum = 0
	for k in range(1, 31):
		k_sum += (1.0/k) * np.sin(k*np.pi*z/H) * \
		np.exp( (-1.0 * (k*k) * (np.pi*np.pi) * v * t / (H*H)))
 
	k_sum *= 2.0 * u_w/np.pi
 
	return u_w * (1.0-(z/H)) - k_sum

def load_mamico_data(csv_file):
	""" Get CSV data from one cycle for mamico
	"""
	# load data in pandas DataFrame
	df = pd.read_csv(csv_file, sep=";")
	df.drop(df.columns[7], axis=1, inplace=True)
	df.columns = ['i','j','k','v_x','v_y','v_z','mass']
	# get Avg x veliocity per z layer
	res = pd.DataFrame(columns=['z', 'v_x', 'mass'],index=index)
	for i,_ in res.iterrows():
		sub_df = df.loc[df.loc[:,'k'] == i]
		res.at[i, 'z'] = 2.5+7.5+1.25+(i-index[0])*2.5
		res.at[i, 'v_x'] = sub_df.loc[:,'v_x'].mean()
		res.at[i, 'mass'] = sub_df.loc[:,'mass'].sum()
	return res

def plot_one_coupling_cycle(cc, color, ax, path):
	""" Plot analytical Couette profile and velocities
	"""
	z = np.linspace(0,H,num=21)
	u = couette_analytic(z, cc*dt)
	ax.plot(z, u, color=color, linewidth=0.5, label='Coupling cycle no. '+str(cc))
	""" Plot MaMiCo velocities
	"""
	csv_file = os.path.join(path, 'CouetteAvgMultiMDCells_0_0_'+str(cc)+'.csv')
	res = load_mamico_data(csv_file)
	ax.plot(res.loc[:,'z'], res.loc[:,'v_x'], marker='.', linestyle='none', markersize=4, color=color)
	""" Mass
	"""
	print('rho:'+str(res.loc[:,'mass'].sum()/volume_inner_md_domain) + ' at ' + str(cc))


def main(argv):
	opts, args = getopt.getopt(argv,'hp::',['path='])
	path = '.'
	for opt, arg in opts:
		if opt == '-h':
			print('profiles.py -p path where path is the folder container MaMiCo outputs')
			sys.exit()
		elif opt in ('-p', '--path'):
			path = arg
	fig = plt.figure(figsize=(7,5))
	ax = fig.add_subplot(111)
	cmap = mpl.colormaps['Spectral']
	norm = mpl.colors.Normalize(vmin=ccs[0], vmax=ccs[-1])

	for cc in ccs:
		plot_one_coupling_cycle(cc, cmap(norm(cc)), ax, path)

	ax.set_xlim(0,H)
	ax.set_xlabel('z', fontsize=14)
	ax.set_ylabel('u_x', fontsize=14)
	ax.set_title(title, fontsize=14)
	ax.grid(visible=True)
	ax.legend()
	plt.show()
	fig.savefig(filename, dpi=250)

if __name__ == "__main__":
   main(sys.argv[1:])
