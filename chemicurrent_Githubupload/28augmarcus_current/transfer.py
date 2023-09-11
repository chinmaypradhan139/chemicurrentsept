import os
import shutil


current=os.getcwd()
dir_name='dU_dG'
os.mkdir(dir_name)
folders=8
main_dir='allinp'

sim_files=['csq_no_','r_no_']
ip_files=['fort.23','fort.45']

for k in range(1,folders+1):
    des=dir_name+'/'+str(k)
    os.mkdir(des)
    dirs=current+'/'+main_dir+'/'+str(k)+'/running'
    runs=current+'/'+main_dir+'/'+str(k)
    for fil in ip_files:
        shutil.copy(runs+'/'+fil,des)
    for sims in sim_files:
        shutil.copy(dirs+'/'+sims+str(k),des)

