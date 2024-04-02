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
import pandas as pd
import datetime as dt
from datetime import timedelta


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

variables = ['dswrf','tasmax', 'tasmin', 'uas', 'vas', 'mrso','cape','pr','tdps','SM_SubX']
var = variables[0]

#Add a new dataset for RZSM and save into home_dir
sub1_dir = f'{subX_dir}/SM_converted_m3_m3'
os.chdir(sub1_dir)

for file in glob('*.nc4'):
    try:
        xr.open_dataset(f'{subX_dir}/RZ{file}')
    except FileNotFoundError:
        
        SubX_file=xr.open_dataset(file)
        
        def date_file_info(SubX_file):
    
            a_date_in= SubX_file.lead.values
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)
    
            return(a_date_out)
    
        julian_list = date_file_info(SubX_file)
        
        #re-assign coords
        SubX_file=SubX_file.assign_coords(lead=julian_list)
        
        #save to home directory
        SubX_file.to_netcdf(f'{subX_dir}/RZ{file}')
