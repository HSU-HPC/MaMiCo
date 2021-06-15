from operator import truediv
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt


if (len(sys.argv) > 1):
    filename=sys.argv[1]
    step=sys.argv[2]
    
else:
    print("Provide filename and timestep")
    exit(1)

actual = pd.read_csv(filename, header=None, sep=';')

actual=actual[actual[0]==int(step)]

res=np.true_divide(np.array(actual[14].values.tolist()), np.array(actual[8].values.tolist()))

towrite=pd.DataFrame(res).transpose()
towrite.to_csv("output.csv", sep=",", header=False, index=False)