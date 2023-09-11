
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#Reading number of cores---------------------------
f=open("fort.23", "r")
i=0

lines=f.readlines()
for line in lines:
    i=i+1
    line=line.strip().split()
    if (i==16):
        cores=int(line[0])
        print(cores)
f.close()
#Averaging-----------------------------------------

init=np.loadtxt(os.path.join("running","1",sys.argv[1]))
data=init[:,1]
for j in range(2,cores+1):
    raw=np.loadtxt(os.path.join("running",str(j),sys.argv[1]))
    data=data+raw[:,1]

data=data/cores
#--------------------------------------------------

plt.plot(init[:,0],1-data)
plt.show()

