import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from pathlib import Path

mpl.use('Agg')
mpl.rcParams['lines.linewidth'] = 2
from matplotlib import pyplot as plt

if (len(sys.argv) > 3):
    csv_file = sys.argv[1]
    csv_type=sys.argv[2]
    stepsize=int(sys.argv[3])
    steps=int(sys.argv[4])
    if not Path(csv_file).is_file():
        print("Not a valid file")
        exit(2)
else:
    print("Provide filename, type (md or macro), step size and numper of plots")
    exit(1)

if(csv_type=="md"):
    lower=4
    upper=10
    z_column=6
    velocity_column=14
    plot_lower=18.75
    plot_upper=23.75+7.5
elif(csv_type=="macro"):
    lower=1
    upper=13
    z_column=3
    velocity_column=8
    plot_lower=50/14*1.5
    plot_upper=50/14*12.5
else:
    print("not a valid type")
    exit(3)




# load data in pandas DataFrame
data = pd.read_csv(csv_file, sep=";", header=None)
if(csv_type=="md"):
    data[8].replace({0 : 1}, inplace=True)
if(csv_type=="macro"):
    data = data[data[1] != 0]
    data = data[data[1] != 13]
    data = data[data[2] != 0]
    data = data[data[2] != 13]
    data = data[data[3] != 0]
    data = data[data[3] != 13]


def couette_analytic(z, t):
    """ Analytic Couette startup equation
    """
    u_w = 1.5
    H = 50.
    v = 2.14/ .813037037

    k_sum = 0
    for k in range(1, 100):
        k_sum += (1.0 / k) * np.sin(k * np.pi * z / H) * \
                 np.exp((-1.0 * (k * k) * (np.pi * np.pi) * v * t / (H * H)))

    k_sum *= 2.0 * u_w / np.pi

    return u_w * (1.0 - (z / H)) - k_sum


def loadAvgDataFromNodeCsv(t):
    """ Get CSV data from one cycle and compute
    the average velocity per layer of cells in z-direction
    """
    
    df=data[data[0]==t]

    # get Avg x velocity per z layer
    avgVelocities = []
    for i in range(lower, upper):
        avg = 0
        mass = 0
        if(csv_type=="md"):
            for _, row in df[df[z_column] == i].iterrows():
                avg += (row[velocity_column])/(row[8])
        else:
            for _, row in df[df[z_column] == i].iterrows():
                avg += (row[velocity_column])
        if df[df[z_column] == i].shape[0] > 0:
            avgVelocities.append(avg / df[df[z_column] == i].shape[0])
    print(avgVelocities)
    return avgVelocities


def plot_one_timestep(t, color, ax):
    """ Plot analytical Couette profile and velocities
    of the inner cells from one coupling cycle
    """
    global simulations
    global offset
    """csv_file = str(simulations) + "/CouetteAvgMultiMDCells_0_" + str(t) + ".csv"
    """
    results = loadAvgDataFromNodeCsv(t)
    #print(results)
    z = np.linspace(0, 50, num=21)
    y = couette_analytic(z, t / 2)
    ax.plot(z, y, "-", color=color)
    ax.plot(np.linspace(plot_lower, plot_upper, num=upper-lower), results, "o", color=color)


def get_cmap(n, name='hsv'):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name."""
    return plt.cm.get_cmap(name, n)

plt.style.use("ggplot")
fig, ax = plt.subplots()

n_plots=5
plot_steps=40
cmap=get_cmap(steps+1)
for i in range(1, steps+1):
    plot_one_timestep(i*stepsize, cmap(i), ax)

#plot_one_timestep(30, "green", ax)
#plot_one_timestep(50, "yellow", ax)
plt.savefig("test.png")


exit(0)
