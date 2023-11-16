#!/usr/bin/env python3
import sys, getopt
import numpy as np
import pandas as pd
import matplotlib as mpl
import os
from matplotlib import pyplot as plt
import seaborn as sns

# simulation parameters
dt = 0.25
ccs = [100, 200, 500]

volume_inner_md_domain = 6*6*6*2.5*2.5*2.5
H = 50
u_w = 0.5
index=[4,5,6,7,8,9]
title = 'Velocity profiles'
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
	# get Avg x velocity per z layer
	res = pd.DataFrame(columns=['z', 'avg_velocities', 'std_velocities', 'mass'],index=index)
	for i,_ in res.iterrows():
		sub_df = df.loc[df.loc[:,'k'] == i]
		res.at[i, 'z'] = sub_df.loc[:,'z'].iat[0]
		res.at[i, 'avg_velocities'] = sub_df.loc[:,'v_x'].mean()
		res.at[i, 'std_velocities'] = sub_df.loc[:,'v_x'].std()
		res.at[i, 'avg_mass'] = sub_df.loc[:,'m'].sum()
	return res

def load_external_data(csv_file):
	""" Get CSV data from one cycle for the external solver
	"""
	# load data in pandas DataFrame
	df = pd.read_csv(csv_file, sep=",")
	return df

def plot_one_coupling_cycle(cc, color, ax, path):
	""" Plot analytical Couette profile and velocities
	"""
	z = np.linspace(0,H,num=21)
	u = couette_analytic(z, cc*dt)
	ax.plot(z, u, color=color, linewidth=0.5, label='Coupling cycle no. '+str(cc))
	""" Plot MaMiCo velocities
	"""
	csv_file = os.path.join(path, 'mamico', 'results_0_'+str(cc)+'.csv')
	res = load_mamico_data(csv_file)
	ax.errorbar(res.loc[:,'z'], res.loc[:,'avg_velocities'], yerr = res.loc[:,'std_velocities'], marker='.', linestyle='none', linewidth=0.5,
    markersize=4, color=color, capsize=0.5, elinewidth=0.5)
	""" Mass
	"""
	#print('z:'+str(res.loc[:,'z']))
	print('rho:'+str(res.loc[:,'avg_mass'].sum()/volume_inner_md_domain))
	""" Plot external velocities
	"""
	csv_file = os.path.join(path, 'openfoam', 'postProcessing', 'sampleDict', str(int(cc*dt)),'l_U.csv')
	res = load_external_data(csv_file)
	ax.plot(res.loc[:,'z'], res.loc[:,'U_0'], color=color, linestyle='none', marker='x', markersize=4);


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
