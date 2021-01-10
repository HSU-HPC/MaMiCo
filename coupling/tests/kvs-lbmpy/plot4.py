#!/usr/bin/env python3

from pandas import read_csv 
import matplotlib.pyplot as mplt

#Variation of helper script plot.py made for comparing exactly four filter .csv outputs.

def plot4(dir, csvs_for_plotting, plot_file_location):
    if len(csvs_for_plotting) != 3:
        print("PLOT4: Incorrect amount of input csvs.")
        return

    mplt.style.use("seaborn")
    fig, ax = mplt.subplots(2,2)
    #unfiltered MD data
    md = read_csv("unfiltered.csv", delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
    ax[0,0].plot(md["Iteration"], md.iloc[:,dir], ".", color="red", label = "MD")

    #all filtered csvs
    #first
    df = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[0].replace("_", " ")
    l = l.replace(".csv", "")
    ax[0,1].plot(df["Iteration"], df.iloc[:,dir], ".", color="green", label = l)
    
    #second
    df = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[1].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,0].plot(df["Iteration"], df.iloc[:,dir], ".", color="orange", label = l)

    #third
    df = read_csv(csvs_for_plotting[2], delimiter=";", usecols=[0,7,8], names=["Iteration", "VelX", "VelY"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[2].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,1].plot(df["Iteration"], df.iloc[:,dir], ".", color="purple", label = l)


    #CS data
    #TODO +3
    cs = read_csv("lbm.csv", delimiter=";", usecols = [0,1,2], names=["Iteration", "VelX", "VelY"])
    for x in range(2):
        for y in range(2):

            #Add CS data to each plot
            ax[x,y].plot(cs["Iteration"], cs.iloc[:,dir], "-", color="blue", label = "Macro solver")

            #get labeling right
            ax[x,y].set_xlabel('coupling cycles')
            ax[x,y].grid(True)
            ax[x,y].legend()
            #TODO: DIRECTION OF VEL
            ax[x,y].set_ylabel('velocity')
    fig.tight_layout()
    #TODO: CUSTOM NAME
    mplt.savefig(plot_file_location)

#YOU HAVE TO CUSTOMIZE THIS
csvs = ["POD.csv", "Gaussian.csv", "NLM.csv"]

#plot for dims 0,1,2 (i.e. x,y,z)
#plot4(0,csvs, "kvs_velx.png")
plot4(1,csvs, "kvs_vely.png")
plot4(2,csvs, "kvs_velz.png")
