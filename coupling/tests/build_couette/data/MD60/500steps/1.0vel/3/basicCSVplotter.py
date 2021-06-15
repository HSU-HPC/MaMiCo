import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

if (len(sys.argv) > 1):
    plotType = sys.argv[1]
    if(plotType=="comparison" or plotType=="unfilteredComparison"):
        if(len(sys.argv)<3):
            print("Provide timestep to compare  to")
            exit(2)
        timestep=int(sys.argv[2])
else:
    print("Provide plotType (unfiltered/averaged/comparison/unfilteredComparison)")
    exit(1)
    
addendum="_clean"

if(addendum=="_clean"):
    index=4
else:
    index=17

df = pd.read_csv('output.csv', header=None)
actual = pd.read_csv('writer1'+addendum+'.csv', header=None, sep=';')
fig, ax = plt.subplots()
#print(df.values.tolist())
if(plotType=="averaged"):
    vlist=df.values.tolist()[0][:216]
    avglist=[]
    print(int(len(vlist)/36))
    for i in range(int(len(vlist)/36)):
        avglist.append(0)
        for j in range(36):
            avglist[i]+=vlist[j+i*36]
        avglist[i]/=36
    print(avglist)
    ax.plot(avglist)
elif(plotType=="unfiltered"):
    ax.plot(df.values.tolist()[0][:216])
elif(plotType=="comparison"):
    vlist=df.values.tolist()[0][:216]
    avglist=[]
    print(int(len(vlist)/36))
    for i in range(int(len(vlist)/36)):
        avglist.append(0)
        for j in range(36):
            avglist[i]+=vlist[j+i*36]
        avglist[i]/=36
    print(avglist)
    ax.plot(avglist)
    
    actualList=actual[actual[0]==timestep].iloc[:,index].tolist()
    actualAvgList=[]
    for i in range(int(len(actualList)/36)):
        actualAvgList.append(0)
        for j in range(36):
            actualAvgList[i]+=actualList[j+i*36]
        actualAvgList[i]/=36
    print(actualAvgList)
    ax.plot(actualAvgList, color="red")
elif(plotType=="unfilteredComparison"):
    ax.plot(df.values.tolist()[0][:216])
    ax.plot(actual[actual[0]==timestep].iloc[:,index].tolist(), color="red")
else:
    print("Provide plotType (unfiltered/averaged/comparison/unfilteredComparison)")

plt.savefig("test.png")
