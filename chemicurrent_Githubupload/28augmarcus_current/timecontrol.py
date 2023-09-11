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
##################################################################
outfile={"fort.103":"M_no_","fort.107":"iM_no_"}
tar="var_dG"
no_fold=50
inputfil="fort.23"
temp_dir="running"
script="./myscript.sh"
folds_in_tar=6
cwd=os.getcwd()
outputfil="MTemp.csv"
##################################################################
def objective(x,a,b):
    return a*x+b

##################################################################

#h=open(outputfil,'w')

#h.write("Gamma"+","+"dG"+","+"Er"+","+"KT"+","+"slope"+","+"chemicurrent (pA)"+","+"% error"+","+"intercept"+"\n") 
#h.close()
for i in range(1,folds_in_tar+1):
    dest=os.path.join(".",tar,str(i))

    shutil.copy("marcus_plot.f90",dest)
    shutil.copy("main_marcus.f90",dest)
    
    f = open(os.path.join(dest,inputfil), "r")
    lines=f.readlines()
    for line in lines:
        line=line.split()
        if line[1]=="!Gamma":
            Gamma=line[0]
        if line[1]=="!dG":
            dG=line[0]
        if line[1]=="!g":
            g=line[0]
        if line[1]=="!Temp":
            KT=line[0]
        if line[1]=="!time_steps":
            nam=line[1]
            stp=line[0] 

    f.close()
   # sys.exit('Done')

    os.chdir(dest)
    os.system("gfortran marcus_plot.f90 main_marcus.f90")
    os.system("./a.out")
    gi=open("fort.24",'r')
    gine=gi.readline()
    value=float(gine)

#    with fileinput.FileInput(inputfil, inplace=True) as file:
#        for line in file:
#            print(line.replace('rep', str(value)), end='')
     
    with fileinput.FileInput(inputfil, inplace=True) as file:
        for line in file:
            print(line.replace(stp+" "+nam, str(value)+" "+nam), end='')

        
    os.chdir(cwd)
    shutil.copy(os.path.join(dest,inputfil),cwd)

    subprocess.call(script)

    for j in range(1,no_fold+1):
        end_path=os.path.join(".",temp_dir,str(j),"ended")
        print(end_path)
        while not os.path.exists(end_path):
            time.sleep(2)

    subprocess.call("pwd")

    for l, (fil,typ) in enumerate(outfile.items()):
        subprocess.call("pwd")
        os.system("python3 avg.py"+" "+str(fil))
        name=typ+str(i)
        os.system("mv output.txt"+ " "+name)
        os.system("mv "+name+ " ./"+tar+"/"+str(i))
    
   # data=np.loadtxt(os.path.join(dest,"iM_no_"+str(i)))
   # nums=data[:,0]
   # ene=data[:,1]
   # popt,pcov=curve_fit(objective,nums,ene)
   # perr=np.sqrt(np.diag(pcov))
   # slope=popt[0]
   # Er=0.5*2000*0.0002**2*float(g)**2
   # hk=open(outputfil,'a')
   # pA=(1.6E-19/2.418E-17)*1.0E+12
   # hk.write(Gamma+","+dG+","+str(Er)+","+KT+","+str(slope)+","+str(slope*pA)+","+str((perr[0]/slope)*100)+","+str(popt[1])+ "\n") 
   # hk.close() 





######################################################################



