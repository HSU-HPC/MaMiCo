import sys
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')
mpl.rcParams['lines.linewidth'] = 2
from matplotlib import pyplot as plt

if (len(sys.argv) > 1):
    simulations = sys.argv[1]
else:
    simulations = "."
offset = 2.5


def couette_analytic(z, t):
    """ Analytic Couette startup equation
    """
    u_w = 1
    H = 50.
    v = 2.14 / .813037037

    k_sum = 0
    for k in range(1, 1001):
        k_sum += (1.0 / k) * np.sin(k * np.pi * z / H) * \
                 np.exp((-1.0 * (k * k) * (np.pi * np.pi) * v * t / (H * H)))

    k_sum *= 2.0 * u_w / np.pi

    return u_w * (1.0 - (z / H)) - k_sum


def loadAvgDataFromNodeCsv(csv_file):
    """ Get CSV data from one cycle and compute
    the average velocity per layer of cells in z-direction
    """
    # load data in pandas DataFrame
    df = pd.read_csv(csv_file, sep=";", header=None)

    # get Avg x velocity per z layer
    avgVelocities = []
    for i in range(4, 10):
        avg = 0
        mass = 0
        for _, row in df[df[2] == i].iterrows():
            avg += (row[3])
        if df[df[2] == i].shape[0] > 0:
            avgVelocities.append(avg / df[df[2] == i].shape[0])
    return avgVelocities


def plot_one_timestep(t, color, ax):
    """ Plot analytical Couette profile and velocities
    of the inner cells from one coupling cycle
    """
    global simulations
    global offset
    csv_file = str(simulations) + "/CouetteAvgMultiMDCells_0_" + str(t) + ".csv"
    results = loadAvgDataFromNodeCsv(csv_file)
    #print(results)
    z = np.linspace(0, 50, num=21)
    y = couette_analytic(z, t / 2)
    ax.plot(z, y, "-", color=color)
    ax.plot(np.linspace(8.75 + offset, 21.25 + offset, num=6), results, "o", color=color)


def get_cmap(n, name='hsv'):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name."""
    return plt.cm.get_cmap(name, n)


plt.style.use("ggplot")
fig, ax = plt.subplots()

n_plots=5
plot_steps=40
cmap=get_cmap(n_plots)
for i in range(1, n_plots):
    plot_one_timestep(i*plot_steps, cmap(i), ax)

#plot_one_timestep(30, "green", ax)
#plot_one_timestep(50, "yellow", ax)
plt.savefig("test.png")


exit(0)
