#!/usr/bin/env python3

from pandas import read_csv 
import matplotlib.pyplot as mplt

def plot(csvs_for_plotting):
    mplt.style.use("seaborn")
    fig, ax = mplt.subplots(2,1)
    for dir in range(2):
        #unfiltered MD data
        df = read_csv("unfiltered.csv", delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
        ax[dir].plot(df["Iteration"], df.iloc[:,dir+1], ".", color="red", label = "MD")

        #all filtered csvs
        for csv in csvs_for_plotting:
            df = read_csv(csv, delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
            #reformat .csv file name
            l = csv.replace("_", " ")
            l = l.replace(".csv", "")
            ax[dir].plot(df["Iteration"], df.iloc[:,dir+1], ".", label = l)

        #CS data
        #TODO +3
        df = read_csv("lbm.csv", delimiter=";", usecols = [0,1,2], names=["Iteration", "VelX", "VelY"])
        ax[dir].plot(df["Iteration"], df.iloc[:,dir+1], "-", color="blue", label = "CS")

        ax[dir].set_xlabel('coupling cycles')
        ax[dir].grid(True)
        ax[dir].legend()
    ax[0].set_ylabel('velocity_x')
    ax[1].set_ylabel('velocity_y')
    fig.tight_layout()
    mplt.savefig("plot.png")

csvs = ["gauss.csv", "median.csv"]
plot(csvs)
