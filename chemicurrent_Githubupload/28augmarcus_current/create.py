import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fileinput
import csv
import sys
import pandas as pd

cols = ["barrier_height","temperature","g","marcus_decay ","% rerror","\n"]

tar="var_dG"
folders=6
inputfil="fort.23"
filtyp=['iM_no_']
outputfil="dG_data.txt"

                                                             
fo = open(outputfil, 'w')
fo.write(",".join(cols))
fo.close()

##################################################################
def objective(x,a,k,b):
    return a*np.exp(-k*x)+b

##################################################################

slope=[]
bar_h=[]
exo=[]
for i in range(1,folders+1):
    for fil in filtyp:
        dest=os.path.join(".",tar,str(i))
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
        fit=[]
        Er=0.5*2000*0.0002**2*float(g)**2
        data=np.loadtxt(os.path.join(dest,fil+str(i)))                        
        nums=data[:,0]                                                         
        ene=data[:,1]                                                          
        guess=[0.02,0.00003,0.99]                                               
        popt,pcov=curve_fit(objective,nums,ene,guess)                          
        err=np.sqrt(np.diag(pcov))                                                                                                                                     
        height=(dG+Er)**2/(4*Er)          
        out=[str(height),KT,g,str(popt[1]),"\n"]
        fk = open(outputfil, 'a')                                                  
        fk.write(",".join(out))                                                    
        fk.close()
        fit=[objective(num,popt[0],popt[1],popt[2]) for num in nums]
        slope.append(popt[1])
        bar_h.append(height)
        exo.append(dG)
 #       plt.plot(nums,fit)
 #       plt.plot(nums,ene)
        #plt.show()
        del fit[:] 


#plt.show()
df=pd.DataFrame({"dG" : exo,
                    "bar_h" : bar_h,
                    "slope" : slope})

df.sort_values(by="bar_h",inplace=True)
print(df)
plt.plot(df["dG"],df["slope"])
plt.show()

