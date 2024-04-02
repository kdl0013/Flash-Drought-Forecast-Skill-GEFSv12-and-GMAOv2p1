#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove unneeded S dimension in SubX files.
Don't run multiprocessing, for some reason it was just messing up files. 
Couldn't figure out why

Also, switch coordinates using nco operators to be on same grid as EMC GEFSv12 and
MERRA 2

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob
import pandas as pd
import datetime as dt
from datetime import timedelta
from multiprocessing import Pool

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='RSMAS'
subX_dir = f'{dir1}/Data/SubX/fromCasper'
out_dir = f'{subX_dir}/S_dim_removed'
os.system(f'mkdir -p {out_dir}')
os.chdir(subX_dir)


#%%
'''Apparently for RSMAS, there are different days of the week that are initialized,
so let us figure out what day of the week each day is initialized'''



dow=[]
for file in sorted(glob('mrso*.nc4')):
    # day_week = file
    dow.append(pd.to_datetime(file.split('.')[0][-10:]).dayofweek)
    
dow[0:100]


'''Yes there are different initialization days for each model.

The solution will be to only keep the file if the S date (within the opened file) is the exact same
as the file date

From manual inspection, the S date is specific to the actual intialization date.

IMPORTANT. The S date is actually file init date. Make sure the file S[1] dim 
is the same as the file name

'''

#First remove the files that aren't good
for file in sorted(glob('*RSMAS*.nc4')):
    open_f = xr.open_dataset(file)
    open_f.close()
    if open_f.S.values[1] == np.datetime64(file.split('.')[0].split('_')[-1]):
        pass
    else:
        os.system(f'rm {file}')

'''All files already have values, this was completed in CASPER cluster and 
The script is in the NCAR_sciprts directory'''

#%%
#Then rename to .nc in accordance with other Subx files
for file in sorted(glob('*RSMAS*.nc4')):
    file_new=file[0:-1]
    os.system(f'mv {file} {file_new}')

