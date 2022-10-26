#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find bad anomaly mean files

@author: kdl
"""
import os
import xarray as xr
from glob import glob
import numpy as np

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
subX_dir = f'{dir1}/Data/SubX/'

model_array=['GMAO', 'ESRL', 'RSMAS', 'ECCC', 'NRL', 'EMC']
var_list=["tasmin" "tasmax" "actual_vapor_pressure" "srad" "windspeed" "tas","pr"]

missing_data = []
for model in model_array:
    print(f'Working on model {model}')
    for var in var_list:
        for file in sorted(glob(f'{subX_dir}/{model}/anomaly/mean_for_ACC/*.nc4')):
            try:
                a=xr.open_dataset(file,engine='netcdf4')
                a.close()
                if len(np.unique(a[list(a.keys())[0]])) < 1000:
                    output = f'{file}'
                    missing_data.append(output)
            except OSError:
                missing_data.append(file)
            except FileNotFoundError:
                missing_data.append(file)
            
'''EMC Penman looks like it may not have data becasue of ncview, 
but it does'''
# a=xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/EMC/anomaly/mean_for_ACC/ETo_mean_11_models_Penman_2000-10-27.nc4')

'''GMAO has issues with precipatiation having not data'''
# a=xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/anomaly/mean_for_ACC/pr_mean_2000-01-16.nc4')
# len(np.unique(a[list(a.keys())[0]]))


#Need to save the missing_data to a txt file
np.savetxt(f'{dir1}/Scripts/missing_data_anomaly.txt', missing_data, fmt='%s'   )
