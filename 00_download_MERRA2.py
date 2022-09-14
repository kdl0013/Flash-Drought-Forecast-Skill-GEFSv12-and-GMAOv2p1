#!/usr/bin/env python3

'''Create wget script for MERRA2 to be run in parallel

How to download:
https://disc.gsfc.nasa.gov/data-access


NOTES:
    This script doesn't download the data as the correct 
    1.) variable name

    
SOLUTION:
    Just merge files after they are downloaded. This will save future storage space.

'''
import os
import datetime as dt
import numpy as np
from glob import glob
import pandas as pd

home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
script_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts'

merra_dir = f'{home_dir}/Data/MERRA2'
os.chdir(f'{merra_dir}')

os.listdir()

download_list = []

file_list_name=["RZSM.txt","albedo_netshortFlux_netdownFlux.txt","temperature.txt"]
count=0
for i_x,f in enumerate(file_list_name):
    
    start_date = pd.to_datetime(dt.date(2000,1,1))
    save_dir = 'all_variables'
    save_dir = f.split(".")[0]
    os.system(f'mkdir -p {save_dir}')
    open_f=np.loadtxt(f, dtype='str')
    
    for i in open_f:
        date_out = f'{start_date.year}-{start_date.month:02}-{start_date.day:02}'
        out_ = f'wget -nc --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies "{i}" -O {merra_dir}/{save_dir}/{save_dir}_{date_out}.nc4 &'

        # os.system(f'{out_}')
        download_list.append(out_)
        start_date = pd.to_datetime(start_date) + np.timedelta64(1,'D')
        count+=1
        if count % 20 ==0:
            download_list.append('wait')

download_list[0]   
   
name = 'wget_MERRA'
np.savetxt(f'{script_dir}/{name}.txt',download_list, fmt="%s")

os.system(f'mv {script_dir}/{name}.txt {script_dir}/{name}.sh')

#save a smaller script. For some reason, it didn't want to save earlier.
name_o = 'wget_MERRA_small'
np.savetxt(f'{script_dir}/{name_o}.txt',download_list[-40:], fmt="%s")

os.system(f'mv {script_dir}/{name_o}.txt {script_dir}/{name_o}.sh')
