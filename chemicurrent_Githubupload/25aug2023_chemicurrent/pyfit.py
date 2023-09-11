import numpy as np
import matplotlib.pyplot as mp
from scipy.optimize import curve_fit

data=np.loadtxt('ir_no_1')

###############################
def objective(x,a,b,k):
    return a*np.exp(-k*x)+b

###############################

time=data[:,0]
pop=data[:,1]

guess=[0.02,0.5,0.00003]
popt,pcov=curve_fit(objective,time,pop,guess)
mp.plot(time,pop)


mp.plot(time,objective(time,popt[0],popt[1],popt[2]))
mp.show()
