#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove unneeded S dimension in SubX files

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob


dir1 = 'main_dir'
mod = 'model'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='GMAO'
subX_dir = f'{dir1}/Data/SubX/{mod}'
os.chdir(subX_dir)

#%%#SMERGE Soil Moisture
#SubX.dat depth file

variables = ['dswrf','tasmax', 'tasmin', 'uas', 'vas', 'mrso','cape','pr','tdps','SM_SubX']
var = variables[0]

for var in variables:
    if var=='SM_SubX':
        sub1_dir = f'{subX_dir}/SM_converted_m3_m3'
        var_name = 'SM_SubX_m3_m3_value'
    else:
        sub1_dir = subX_dir
        var_name=var
    os.chdir(sub1_dir)
    try:
        xr.open_dataset(sorted(glob(f'{var}*.nc4'))[0])
    except IndexError:
        for file in sorted(glob(f'{var}*.nc')):
            open_f = xr.open_dataset(file)
            open_f[f'{var_name}'][1,:,:,:,:] = np.nan
            #Drop S dimension to save storage space
            open_f = open_f.dropna(dim='S',how='all')
            open_f.close()
            open_f.to_netcdf(f'{file}4')
            open_f.close()
            os.system(f'rm {file}')
        