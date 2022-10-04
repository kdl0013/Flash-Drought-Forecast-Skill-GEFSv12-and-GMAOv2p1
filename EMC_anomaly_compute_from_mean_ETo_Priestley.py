#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly for ETo (no RZSM). Just subtract mean from file.
    




@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn

# dir1 = 'main_dir'
# start_ = int('start_init')
# end_ = start_ + int('init_step')
# model_NAM1 = 'EMC'
# var = 'variable_'  #for anomaly function

# Test for 1 step size and model
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'EMC'
name_ = 'Priestley'
#%%
var = f'ETo_{name_}'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
anom_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
mean_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly/mean_for_ACC'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
mean_ensemble_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly/MEM'
os.system(f'mkdir -p {mean_ensemble_dir}')

script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)
#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{home_dir}/{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var)    

'''Steps for anomaly calculation. 
Subtract mean from the file
        
'''
print('Creating anomaly files. Saving into multi ensemble mean (MEM) and as individuals models')

for file in sorted(glob(f'{var}*.nc4')):
    open_f = xr.open_dataset(file)
    open_f.close()
    
    if model_NAM1 == 'EMC':
        mean_date = f"ETo_mean_11_models_{name_}_2000-{file.split('_')[-1].split('.')[0][5:]}.nc4"
    else:
        mean_date = f"ETo_mean_{name_}_2000-{file.split('_')[-1].split('.')[0][5:]}.nc4"
    open_mean = xr.open_dataset(f'{mean_dir}/{mean_date}')
    open_mean.close()
    #for some reason, you have to convert the datasets. Can't just subtract
    if open_f.M.shape[0] != open_mean.M.shape[0]:
        min_ = min(open_f.M.shape[0],open_mean.M.shape[0])
        anomaly = (open_f.ETo[:,0:min_,:,:,:].to_numpy() - open_mean.ETo_mean[:,0:min_,:,:,:].to_numpy())
    else:
        anomaly = (open_f.ETo[:,:,:,:,:].to_numpy() - open_mean.ETo_mean[:,:,:,:,:].to_numpy())
        min_ = open_f.M.shape[0]
        
    anomaly_MEM = (open_f.ETo[:,0:min_,:,:,:].mean(axis=1).to_numpy() - open_mean.ETo_mean[:,0:min_,:,:,:].mean(axis=1).to_numpy())
    
    if open_f.M.shape[0] < open_mean.M.shape[0]:
        open_f = xr.zeros_like(open_f.isel(M=slice(0,open_f.M.shape[0])))
    elif open_mean.M.shape[0] < open_f.M.shape[0]:
        open_f = xr.zeros_like(open_f.isel(M=slice(0,open_mean.M.shape[0])))

        
    open_f.ETo[:,:,:,:,:] = anomaly
    open_f = open_f.rename(ETo='ETo_anom')
    open_f = open_f.assign_coords(L=np.arange(open_f.L.shape[0])) #now save as lead date for anomaly correlation
    open_f.to_netcdf(f"{anom_dir}/ETo_{name_}_anomaly_{model_NAM1}_{file.split('_')[-1]}")
    os.system(f"ncks -4 -L 4 {anom_dir}/ETo_{name_}_anomaly_{model_NAM1}_{file.split('_')[-1]} {anom_dir}/a_ETo_{name_}_anomaly_{model_NAM1}_{file.split('_')[-1]}")
    os.system(f"mv {anom_dir}/a_ETo_{name_}_anomaly_{model_NAM1}_{file.split('_')[-1]} {anom_dir}/ETo_{name_}_anomaly_{model_NAM1}_{file.split('_')[-1]}")
    #make new dataset for anomaly multi_ensemble mean (MEM)
    
    #Convert to an xarray object
    var_OUT = xr.Dataset(
        data_vars = dict(
            ETo_anom_MEM = (['S','L','Y','X'],  anomaly_MEM),
        ),
        coords = dict(
            X = open_f.X.values,
            Y = open_f.Y.values,
            L = open_f.L.values,
            S = np.atleast_1d(pd.to_datetime(file.split('_')[-1].split('.')[0]))
        ),
        attrs = dict(
            Description = 'Reference crop evapotranspiration (mm/day). Priestley-Taylor formula'),
    )              
    var_OUT.to_netcdf(f"{mean_ensemble_dir}/ETo_{name_}_anomaly_MEM_{model_NAM1}_{file.split('_')[-1]}")
    os.system(f"ncks -4 -L 4 {mean_ensemble_dir}/ETo_{name_}_anomaly_MEM_{model_NAM1}_{file.split('_')[-1]} {mean_ensemble_dir}/a_ETo_{name_}_anomaly_MEM_{model_NAM1}_{file.split('_')[-1]} ")
    os.system(f"mv {mean_ensemble_dir}/a_ETo_{name_}_anomaly_MEM_{model_NAM1}_{file.split('_')[-1]}  {mean_ensemble_dir}/ETo_{name_}_anomaly_MEM_{model_NAM1}_{file.split('_')[-1]}")

