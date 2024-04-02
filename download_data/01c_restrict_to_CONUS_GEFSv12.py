#!/usr/bin/env python3

'''
Confine to CONUS region.

'''
import os
import datetime as dt
import numpy as np
import xarray as xr
from glob import glob
import pandas as pd


dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'

# #for linux
in_dir = f'{dir1}/Data/SubX/EMC/'
os.chdir(in_dir)

save_out_dir_GEFS = f'{dir1}/Data/SubX/EMC/CONUS'
os.system(f'mkdir {save_out_dir_GEFS}')


all_files = os.listdir()

for file in all_files:
    try:
        xr.open_dataset(f'{save_out_dir_GEFS}/{file}')
    except FileNotFoundError:
        os.system(f'ncks -d X,235.0,293.0 -d Y,24.0,50.0 {file} {save_out_dir_GEFS}/{file}')
        os.system(f'ncks -4 -L 1 {save_out_dir_GEFS}/{file} {save_out_dir_GEFS}/{file}4')
        os.system(f'rm {save_out_dir_GEFS}/{file}')


# for i in vars_to_process:
#     os.system(f'mkdir -p {out_dir}/{i}')
#%%

