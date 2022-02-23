import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['lines.linewidth'] = 2
from matplotlib import pyplot as plt
 

offset = 2.5
 
def couette_analytic(z, t):
	""" Analytic Couette startup equation
	"""
	u_w = 1.5
	H = 50.
	v = 2.14 / .813037037
 
	k_sum = 0
	for k in range(1, 1001):
		k_sum += (1.0/k) * np.sin(k*np.pi*z/H) * \
		np.exp( (-1.0 * (k*k) * (np.pi*np.pi) * v * t / (H*H)))
 
	k_sum *= 2.0 * u_w/np.pi
 
	return u_w * (1.0-(z/H)) - k_sum

def loadAvgDataFromNodeCsv(csv_file):
	""" Get CSV data from one cycle and compute
	the average velocity per layer of cells in z-direction
	"""
	# load data in pandas DataFrame
	df = pd.read_csv(csv_file, sep=";", header=None)
 
	# get Avg x velocity per z layer
	avgVelocities = []
	for i in range(4,10):
		avg = 0
		mass = 0
		j = 0
		for _,row in df[df[2] == i].iterrows():
			avg += (row[3])
			#break
		if df[df[2] == i].shape[0] > 0:
			avgVelocities.append(avg/df[df[2] == i].shape[0])
		#avgVelocities.append(avg)
	return avgVelocities
 
def plot_one_timestep(t, color, ax):
	""" Plot analytical Couette profile and velocities
	of the inner cells from one coupling cycle
	"""
	global simulations
	global offset
	csv_file = "CouetteAvgMultiMDCells_0_" + str(t) + ".csv"
	results = loadAvgDataFromNodeCsv(csv_file)
	print(results)
	z = np.linspace(0,50,num=21)
	y = couette_analytic(z, t/4) #multiply MD timestep by number of MD per coupling, multiply that factor with t here
	ax.plot(z, y, "-", color=color)
	ax.plot(np.linspace(7.5+1.25+offset, 7.5+1.25+12.5+offset, num=6), results, "o", color=color)
 
plt.style.use("ggplot")
fig, ax = plt.subplots()
	 
plot_one_timestep(30, "green", ax)
plot_one_timestep(50, "yellow", ax)
plot_one_timestep(80, "blue", ax)
plot_one_timestep(130, "cyan", ax)
plot_one_timestep(250, "orange", ax)
plot_one_timestep(500, "darkred", ax)
plt.savefig("test.png", format="png")
plt.show()
 
exit(0)
