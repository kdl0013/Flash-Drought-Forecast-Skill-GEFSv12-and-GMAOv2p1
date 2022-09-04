#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert GMAO soil moisture from kg/m2 to m3/m3. We are given depth. Just take
kg/m2 / m * 1000

@author: kdl
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
from datetime import timedelta
import pandas as pd
from glob import glob
from multiprocessing import Pool


# dir1 = 'main_dir'
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

model='GMAO'
subX_dir = f'{dir1}/Data/SubX/{model}'

os.chdir(f'{subX_dir}')

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'mrso')   

for file in sorted(glob('mrso*')):
    try:
        open_f = xr.open_dataset(file)
        
        file_date = np.datetime64(file.split('.')[0].split('_')[-1])
        
        if open_f[list(open_f.keys())[0]].S.values[0] != file_date:
            os.system(f'rm *{file_date}*')
    except ValueError:
        pass #no data in file

#%%
