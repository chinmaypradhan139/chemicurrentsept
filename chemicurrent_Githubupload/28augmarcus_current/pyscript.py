
import os
import shutil
import time
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fileinput


subprocess.call("(rm output.txt)",shell=True)
subprocess.call("(python3 avg.py fort.107)",shell=True)


f = open('fort.23', "r")
lines=f.readlines()
for line in lines:
    line=line.split()
    if line[1]=="!Gamma":
        Gamma=line[0]
    if line[1]=="!dG":
        dG=float(line[0])
    if line[1]=="!g":
        g=float(line[0])
    if line[1]=="!Temp":
        KT=line[0]
    if line[1]=="!time_steps":
        nam=line[1]
        stp=line[0] 

f.close()

Er=0.5*2000*(0.0002)**2*g**2

height=(dG+Er)**2/(4*Er)


##################################################################
def objective(x,a,b,k):
    return a*np.exp(-k*x)+b

##################################################################
data=np.loadtxt("output.txt")
nums=data[:,0]
ene=data[:,1]
guess=[0.02,0.5,0.00003]
popt,pcov=curve_fit(objective,nums,ene,guess)
cols = ["Gamma","g","dG","barrier_height","KT","rate constant","error","population difference","\n"]
error=np.sqrt(np.diag(pcov))

if 
    fo = open('para.txt', 'a')
    fo.write(",".join(cols))
    fo.close()



out=[Gamma,str(g),str(dG),str(height),KT,str(popt[2]),str(error[0]),str(ene[-1]-ene[0]),"\n"]


fk = open('para.txt', 'a')
fk.write(",".join(out))
fk.close()

