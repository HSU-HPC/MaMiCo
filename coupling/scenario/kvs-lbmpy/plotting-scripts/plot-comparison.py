#!/usr/bin/env python3

import math

import matplotlib.pyplot as mplt
from pandas import read_csv

#Variation of helper script plot4.py for cases in which CS data is not present
#This implies no SNR computation

def plot_comparison(dir, csvs_for_plotting, plot_file_location):
    if len(csvs_for_plotting) != 2:
        print("PLOT-COMPARISON: Incorrect amount of input csvs.")
        return

    mplt.style.use("seaborn")
    fig, ax = mplt.subplots(2,2)

    #unfiltered MD data
    p_md = read_csv(unfiltered_csv, delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_md = read_csv(unfiltered_csv, delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    ax[0,0].plot(p_md["Iteration"], (p_md.iloc[:,dir]) / (m_md.iloc[:,1]), ".", color="red", label = "MD")


    #first filtered dataset
    p_fst = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_fst = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l1 = csvs_for_plotting[0].replace("_", " ")
    l1 = l1.replace("../", "")
    l1 = l1.replace(".csv", "")
    ax[0,1].plot(p_fst["Iteration"], (p_fst.iloc[:,dir]) / (m_fst.iloc[:,1]), ".", color="green", label = l1)

    
    #second filtered dataset
    p_snd = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_snd = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l2 = csvs_for_plotting[1].replace("_", " ")
    l2 = l2.replace("../", "")
    l2 = l2.replace(".csv", "")
    ax[1,0].plot(p_snd["Iteration"], (p_snd.iloc[:,dir]) / (m_snd.iloc[:,1]), ".", color="orange", label = l2)


    #both at the same time
    ax[1,1].plot(p_fst["Iteration"], (p_fst.iloc[:,dir]) / (m_fst.iloc[:,1]), ".", color="green", label = l1)
    ax[1,1].plot(p_snd["Iteration"], (p_snd.iloc[:,dir]) / (m_snd.iloc[:,1]), ".", color="orange", label = l2)

    for x in range(2):
        for y in range(2):

            #get labeling right
            ax[x,y].set_xlabel('coupling cycles')
            ax[x,y].grid(True)
            ax[x,y].legend()
            
            #too lazy for a case-switch
            if dir == 1:
                ax[x,y].set_ylabel('X velocity')
            if dir == 2:
                ax[x,y].set_ylabel('Y velocity')
            if dir == 3:
                ax[x,y].set_ylabel('Z velocity')

    fig.tight_layout()
    mplt.savefig(plot_file_location)

#YOU HAVE TO CUSTOMIZE THIS
unfiltered_csv = "../unfiltered.csv"
csvs = ["../gaussXYZ-python.csv", "../gaussXYZ-cpp.csv"]

#plot for dims 1,2,3 (i.e. x,y,z)
plot_comparison(1,csvs, "kvs_velx_gausscomp_mirror_XYZ.png")
plot_comparison(2,csvs, "kvs_vely_gausscomp_mirror_XYZ.png")
