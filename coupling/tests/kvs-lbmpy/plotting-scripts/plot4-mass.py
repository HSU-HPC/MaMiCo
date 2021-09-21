#!/usr/bin/env python3

from pandas import read_csv 
import math
import matplotlib.pyplot as mplt

#Variation of helper script plot.py made for comparing exactly four filter .csv outputs.
#This version of plot4 works on mass only. It does not compute an SNR or use CS data.

def plot4Mass(csvs_for_plotting, plot_file_location):
    if len(csvs_for_plotting) != 3:
        print("PLOT4: Incorrect amount of input csvs.")
        return

    mplt.style.use("seaborn")
    fig, ax = mplt.subplots(2,2)
    #unfiltered MD data
    m_md = read_csv(unfiltered_csv, delimiter=";", usecols=[0,7], names=["Iteration", "Mass"], index_col=None)
    ax[0,0].plot(m_md["Iteration"], m_md.iloc[:,1], ".", color="red", label = "MD")
    

    #all filtered csvs
    #first
    m_fst = read_csv(csvs_for_plotting[0], delimiter=";", usecols=[0,7], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[0].replace("_", " ")
    l = l.replace(".csv", "")
    ax[0,1].plot(m_fst["Iteration"], m_fst.iloc[:,1], ".", color="green", label = l)
    
    
    #second
    m_snd = read_csv(csvs_for_plotting[1], delimiter=";", usecols=[0,7], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[1].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,0].plot(m_snd["Iteration"], m_snd.iloc[:,1], ".", color="orange", label = l)
    

    #third
    m_thrd = read_csv(csvs_for_plotting[2], delimiter=";", usecols=[0,7], names=["Iteration", "Mass"], index_col=None)
    #reformat .csv file name
    l = csvs_for_plotting[2].replace("_", " ")
    l = l.replace(".csv", "")
    ax[1,1].plot(m_thrd["Iteration"], m_thrd.iloc[:,1], ".", color="purple", label = l)
    

    for x in range(2):
        for y in range(2):
      
            #get labeling right
            ax[x,y].set_xlabel('coupling cycles')
            ax[x,y].grid(True)
            ax[x,y].legend()

            ax[x,y].set_ylabel('Mass')
            ax[x,y].set_ylim(-0.5,1.5)

    fig.tight_layout()
    mplt.savefig(plot_file_location)

#YOU HAVE TO CUSTOMIZE THIS
unfiltered_csv = "../unfiltered.csv"
csvs = ["../POD.csv", "../Gaussian.csv", "../NLM_junction.csv"]

#plot for dims 1,2,3 (i.e. x,y,z)
plot4Mass(csvs, "kvs_mass.png")
