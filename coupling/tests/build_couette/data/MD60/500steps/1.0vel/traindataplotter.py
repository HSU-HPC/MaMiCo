import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

def get_cmap(n, name='hsv'):
    """Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name."""
    return plt.cm.get_cmap(name, n)

n=4
timestep=250

for gauss in range(1,5):
    actual=[]
    for i in range(1,n+1):
        actual.append(pd.read_csv(str(i)+'/writer60_after'+str(gauss)+'.csv', header=None, sep=';'))
    
    fig, ax = plt.subplots()
    
    cmap=get_cmap(n)
    
    for i in range(n):
        print(np.true_divide(np.array(actual[i][actual[i][0]==timestep].iloc[:,14].values.tolist()), np.array(actual[i][actual[i][0]==timestep].iloc[:,8].values.tolist())))
        ax.plot(np.true_divide(np.array(actual[i][actual[i][0]==timestep].iloc[:,14].values.tolist()), np.array(actual[i][actual[i][0]==timestep].iloc[:,8].values.tolist())), color=cmap(i))
    
    plt.savefig("after"+str(gauss)+".png")









