from operator import truediv
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt


actual = pd.read_csv('writer2.csv', header=None, sep=';')
actual=actual[actual[1]==3]
actual=actual[actual[2]==3]
actual=actual[actual[3]!=0]
actual=actual[actual[3]!=1]
actual=actual[actual[3]!=2]
actual=actual[actual[3]!=3]

actual=actual[actual[3]!=10]
actual=actual[actual[3]!=11]
actual=actual[actual[3]!=12]
actual=actual[actual[3]!=13]

actual=actual[[0,3,8]]

actual=actual.reset_index(drop=True)
#relevant=actual[actual[0]==1].iloc[:,17].values.tolist()

#for i in range(1, actual.iloc[[-1]][0])

newarray= np.empty((0, 5), float)
end=int(actual.iloc[[-1]][0])+1

for step in range(1, end):
    for i in range(4,10):
        value=float(actual[(actual[0]==step) & (actual[3]==i)].iloc[[0]][8])
        for j in range(4, 10):
            for k in range(4, 10):
                newarray = np.append(newarray, np.array([[step, k, j, i, value]]), axis=0)

towrite=pd.DataFrame(newarray)
towrite = towrite.astype({0: int,1:int,2:int,3:int})
towrite.to_csv("CleanTraining250.csv", sep=";", header=False, index=False)


