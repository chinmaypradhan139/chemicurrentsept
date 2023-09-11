import os
import shutil
import time
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fileinput
import csv
import math as mp
##################################################################
decor=int(sys.argv[1])
outfile1={"fort.101":"csq_no_","fort.102":"icsq_no_","fort.104":"r_no_","fort.106":"ir_no_"}
outfile2={"fort.101":"dcsq_no_","fort.102":"idcsq_no_","fort.104":"dr_no_","fort.106":"idr_no_"}
tar="temp_var"

inputfil="fort.23"
temp_dir="running"
script="./myscript.sh"
folds_in_tar=6
cwd=os.getcwd()
print(decor,"decor")
if decor==2:
    outputfil="var_temp.csv"
else:
    outputfil="var_temp.csv"
##################################################################
def objective(x,a,b,k):
    return a*np.exp(-k*x)+b

##################################################################

print(outputfil, "outputfil")
cols = ["barrier_height","rate constant mod c**2","% ierror","longtime pop","rate constant traj based ","% rerror","longtime pop","\n"]

fo = open(outputfil, 'w')
fo.write(",".join(cols))
fo.close()


for i in range(1,folds_in_tar+1):
    dest=os.path.join(".",tar,str(i))

    shutil.copy("marcus_plot.f90",dest)
    shutil.copy("main_marcus.f90",dest)
    
    print("1 hi")
    os.system("cat fort.23")

    os.chdir(dest)
    os.system("pwd")
    os.system("gfortran marcus_plot.f90 main_marcus.f90")
    os.system("./a.out")
    gi=open("runtime",'r')
    gine=gi.readline()
    value=float(gine)
    print(value, "value") 

    os.chdir(cwd)
    f = open(os.path.join(dest,inputfil), "r")
    lines=f.readlines()
    for line in lines:
        line=line.split()
        if line[1]=="!Gamma":
            Gamma=line[0]
        if line[1]=="!dG":
            dG=float(line[0])
        if line[1]=="!g":
            g=line[0]
        if line[1]=="!Temp":
            KT=line[0]
        if line[1]=="!time_steps":
            nam=line[1]
            stp=line[0]
        if line[1]=="!1":
            namd=line[1]
            dec=line[0]
        if line[1]=="!num_cores":
            no_fold=int(line[0])
            
    f.close()

    Er=0.5*2000*0.0002**2*float(g)**2

    shutil.copy(os.path.join(dest,inputfil),cwd)

    with fileinput.FileInput(inputfil, inplace=True) as file:
        for line in file:
            print(line.replace(stp+" "+nam, str(value)+" "+nam), end='')
            
    with fileinput.FileInput(inputfil, inplace=True) as file:
        for line in file:
            print(line.replace(dec+" "+namd, str(decor)+" "+namd), end='')
 

    print("4 hi")
    os.system("cat fort.23")
    subprocess.call(script)
    
    for j in range(1,no_fold+1):
        end_path=os.path.join(".",temp_dir,str(j),"ended")
        print(end_path)
        while not os.path.exists(end_path):
            time.sleep(1)
    
    print("hi")
    if decor==2:
        outfile=outfile1
    else:
        outfile=outfile2

    for l, (fil,typ) in enumerate(outfile.items()):
        os.system("python3 avg.py"+" "+str(fil))
        name=typ+str(i)
        os.system("mv output.txt"+ " "+name)
        os.system("mv "+name+ " ./"+tar+"/"+str(i))
    
    
#    if decor==2:
#        ifile=["icsq_no_","ir_no_"]
#    else:
#        ifile=["idcsq_no_","idr_no_"]

#    rate=[]
#    perr=[]
#    inter=[]
#    for ifil in ifile:
#        data=np.loadtxt(os.path.join(dest,ifil+str(i)))
#        nums=data[:,0]
#        ene=data[:,1]
#        guess=[0.02,0.5,0.00003]
#        popt,pcov=curve_fit(objective,nums,ene,guess)
#        err=np.sqrt(np.diag(pcov))
#        rate.append(popt[2])
#        perr.append(err[2])
#        inter.append(popt[1])

#    height=(dG+Er)**2/(4*Er)
#    print(type(height))
#    out=[str(height),str(rate[0]),str(perr[0]),str(inter[0]),str(rate[1]),str(perr[1]),str(inter[1]),"\n"]
#    fk = open(outputfil, 'a')
#    fk.write(",".join(out))
#    fk.close()






######################################################################



