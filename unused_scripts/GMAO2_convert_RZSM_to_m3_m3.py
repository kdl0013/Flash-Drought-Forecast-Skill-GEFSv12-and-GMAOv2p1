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
subX_SM_out_dir = f'{subX_dir}/SM_converted_to_m3_m3'

os.system(f'mkdir -p {subX_SM_out_dir}')

os.chdir(f'{subX_dir}')

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'mrso')   
init_date_list[0] 
init_date_list[-1]

#%%#SMERGE Soil Moisture
#SubX depth file (already converted)
global depth_file
depth_file = xr.open_dataset(f'{dir1}/Data/SubX/SubX_sm_grid_conversion/GMAO_1x1_grid.nc')

def convert_SubX_to_m3_m3(_date):
# for _date in init_date_list:    
    # 
    try:
       xr.open_dataset(f'{subX_SM_out_dir}/SM_{model}_m3_m3_{_date}.nc4')
       print(f'Already completed date {_date}. Saved in {subX_SM_out_dir}')
       
    except FileNotFoundError:
        
        print(f'Working on initialized day {_date} to convert SubX soil moisture into m3/m3 and saving into {subX_SM_out_dir}.')

        SubX_file = xr.open_dataset(f'{subX_dir}/mrso_GMAO_{_date}.nc4')
        var='RZSM'
        
        def date_file_info(SubX_file):

            a_date_in= SubX_file.L.values
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

            return(a_date_out)
        
        julian_list = date_file_info(SubX_file)
        
        out_sm_SubX = xr.zeros_like(SubX_file)
        
        for mod in range(SubX_file.M.shape[0]):
            for lead in range(SubX_file.L.shape[0]):
                # print(mod,lead)
                out_sm_SubX[list(out_sm_SubX.keys())[0]][0,mod, lead, :, :] = SubX_file[list(SubX_file.keys())[0]][0,mod, lead, :, :].values / ((depth_file.gmao_depth[:,:].values)*1000)

        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
        data_vars = dict(
            SM_SubX_m3_m3_value = (['S','model','lead','Y','X'], out_sm_SubX[list(out_sm_SubX.keys())[0]].values),
        ),
        coords = dict(
        S = out_sm_SubX.S.values,
        X = out_sm_SubX.X.values,
        Y = out_sm_SubX.Y.values,
        lead = julian_list,
        model = out_sm_SubX.M.values,
        ),
        attrs = dict(
        Description = f'SubX SM (kg/m2) converted to (m3/m3) based on profile depth {model}.'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{subX_SM_out_dir}/SM_{model}_m3_m3_{_date}.nc4', mode ='w')
        print(f'Completed SM_m3_m3_{_date}.')


#%%
if __name__ == '__main__':
    p = Pool(8)
    p.map(convert_SubX_to_m3_m3,init_date_list)