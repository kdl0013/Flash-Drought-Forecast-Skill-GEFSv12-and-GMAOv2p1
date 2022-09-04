#!/usr/bin/env python3

'''
Some GEFSv12 years (after 2020), have 2 dimensions.
Delete them

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
save_out_dir_GEFS = f'{dir1}/Data/SubX/EMC/'
os.chdir(save_out_dir_GEFS)

vars_to_process= ['pr','soilw1', 'soilw2','soilw3', 'soilw4', 'soilw_bgrnd', 'dlwrf','dswrf','tas']


# for i in vars_to_process:
#     os.system(f'mkdir -p {out_dir}/{i}')
#%%

#Open file, find out if it has 2 values. If it does, take the 2nd value

for var in vars_to_process:
    # var=vars_to_process[0]
    print(f'Working on variable {var} to remove S dimension.')
    all_files = sorted(glob(f'{var}*.nc'))
    for file in all_files:
        open_f = xr.open_dataset(file)
        varname=list(open_f.keys())[0]
        if len(open_f[f'{varname}'].S.values) == 2:
            open_f = open_f.dropna(dim='S',how='all')
            open_f.close()
            len(np.unique(open_f[f'{varname}'].values))
            open_f.to_netcdf(f'{file}4')
            os.system(f'rm {file}')
            os.system(f'mv {file}4 {file}')

