import numpy as np
import matplotlib.pyplot as plt
import os


# get number of metal orbitals----------------------
f=open("fort.23", "r")
i=0

lines=f.readlines()
for line in lines:
    i=i+1
    line=line.strip().split()
    if (i==1):
        orbs=int(line[0])-1
    if (i==16):
        cores=int(line[0])
f.close()
#Averaging-----------------------------------------

init=np.loadtxt(os.path.join("running","1","fort.432"))
data=init[:,1]
for j in range(2,cores+1):
    raw=np.loadtxt(os.path.join("running",str(j),"fort.432"))
    data=data+raw[:,1]

data=data/cores
#get fermi_distribution data--------------------------------------------

#fdist=np.loadtxt(os.path.join("running","1","fort.67"))
#plt.plot(fdist[:,0],fdist[:,1])
#plt.savefig("snaps/1/"+"1"+".png")


#-------------------------------------------------
snaps=int(len(data)/orbs)

no_of_snaps=range(1,snaps+1)

snap_orbs=[y*orbs for y in no_of_snaps]

snap_orbs.insert(0,0)

ene=init[:,0]

energy=ene[0:orbs]
print(snap_orbs)

for j in range(len(snap_orbs)-1):
#    print(j)
    time_snap=data[snap_orbs[j]:snap_orbs[j+1]]
    plt.plot(energy,time_snap)
    plt.savefig("snaps/1/"+str(j)+".png")



