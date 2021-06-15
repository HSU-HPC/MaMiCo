from operator import truediv
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt


if (len(sys.argv) > 3):
    file=sys.argv[1]
    plotType = sys.argv[2]
    timestep=int(sys.argv[3])
else:
    print("Provide filename, plotType (moment/vel) and timestep")
    exit(1)


actual = pd.read_csv(file, header=None, sep=';')
fig, ax = plt.subplots()
#print(df.values.tolist())
if(plotType=="vel"):
    ax.plot(actual[actual[0]==timestep].iloc[:,17].tolist(), color="red")
elif(plotType=="moment"):
    ax.plot(np.array(actual[actual[0]==timestep].iloc[:,14].tolist())/np.array(actual[actual[0]==timestep].iloc[:,8].tolist()), color="red")
else:
    print("unknown type")

plt.savefig(file+".png")
