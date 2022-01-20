#!/usr/bin/env python3

from pandas import read_csv 
import math
import matplotlib.pyplot as mplt

#Variation of helper script plot4.py for cases in which CS data is not present
#This implies no SNR computation

def plot4(dir, csvs_for_plotting, plot_file_location):
    if len(csvs_for_plotting) != 3:
        print("PLOT4: Incorrect amount of input csvs.")
        return

    mplt.style.use("seaborn")
    fig, ax = mplt.subplots(2,2)
    #unfiltered MD data
    p_md = read_csv(unfiltered_csv, delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_md = read_csv(unfiltered_csv, delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    ax[0,0].plot(p_md["Iteration"], (p_md.iloc[:,dir]) / (m_md.iloc[:,1]), ".", color="red", label = "MD")


    #all filtered csvs
    #first
    p_fst = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_fst = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[0].replace("_", " ")
    l = l.replace(".csv", "")
    ax[0,1].plot(p_fst["Iteration"], (p_fst.iloc[:,dir]) / (m_fst.iloc[:,1]), ".", color="green", label = l)

    
    #second
    p_snd = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_snd = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[1].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,0].plot(p_snd["Iteration"], (p_snd.iloc[:,dir]) / (m_snd.iloc[:,1]), ".", color="orange", label = l)


    #third
    p_thrd = read_csv(csvs_for_plotting[2], delimiter=";", usecols=[0,5,6,7], names=["Iteration", "MomentumX", "MomentumY", "MomentumZ"], index_col=None)
    m_thrd = read_csv(csvs_for_plotting[2], delimiter=";", usecols=[0,4], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[2].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,1].plot(p_thrd["Iteration"], (p_thrd.iloc[:,dir]) / (m_thrd.iloc[:,1]), ".", color="purple", label = l)

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
csvs = ["../unfiltered.csv", "../gaussXYZ-python.csv", "../gaussXYZ-cpp.csv"]

#plot for dims 1,2,3 (i.e. x,y,z)
plot4(1,csvs, "kvs_velx.png")
plot4(2,csvs, "kvs_vely.png")
#plot4(3,csvs, "kvs_velz.png")
