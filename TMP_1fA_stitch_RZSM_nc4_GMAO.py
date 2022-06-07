#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate reference evapotranspiration (ETo) using ASCE Penman Monteith
    method https://pypi.org/project/refet/ python package : refet
    
    gridMET calculate reference evapotranspiration with the same python package
    
    Calculate EDDI from SubX data (EDDI is already downloaded for historical data)
    


@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod = 'GMAO'
var = 'RZSM'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#create dictionary of with 4 models (GMAO specifically)
model_dirs = {}

if var == 'EDDI':
    for mod in [0,1,2,3]:
        model_dirs[f'Model {mod}'] = f'{home_dir}/EDDI_mod{mod}'
else:
    for mod in [0,1,2,3]:
        model_dirs[f'Model {mod}'] = f'{home_dir}/{var}_anomaly_mod{mod}'

# model_0 = f'{home_dir}/EDDI_mod0/already_completed'
# model_1 = f'{home_dir}/EDDI_mod1/already_completed'
# model_2 = f'{home_dir}/EDDI_mod2/already_completed'
# model_3 = f'{home_dir}/EDDI_mod3/already_completed'


###Files
file_list = os.listdir()

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
init_date_list = return_date_list(var = 'ETo')    

_date=init_date_list[0]
print(f'Saving {var} netcdf empty files into multiple directories')
for _date in init_date_list:
    
    if var == 'EDDI':
        filename = 'EDDI'
        desc = 'Evaporative Demand Drought Index'
    elif var == 'RZSM':
        filename = f'{var}_anomaly'
        desc = 'RZSM anomaly. 7-day average mean by julian day.'
    elif var == 'ETo':
        filename = f'{var}_anomaly'
        desc = 'ETo anomaly. 7-day average mean by julian day.'
    
    template = xr.open_dataset(f'ETo_{_date}.nc4')

    #initialize empty file
    var_OUT = xr.zeros_like(template)
    
    #Add data from each model into the empty file
    for mod in model_dirs.items():
        var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = np.load(f'{filename}_{_date}.npy',allow_pickle=True)
        #print(mod)
        
    #save file as EDDI as netcdf
    var_final = xr.Dataset(
        data_vars = dict(
            Variables = (['S','lead','Y','X'], var_OUT.ETo[:,0,:,:,:].values),
        ),
        coords = dict(
            X = var_OUT.X.values,
            Y = var_OUT.Y.values,
            lead = var_OUT.lead.values,
            S = template.S.values
        ),
        attrs = dict(
            Description = f'{desc}'),
    )                    
    
    #Save as a netcdf for later processing
    
    for mod in model_dirs.items():
        var_final.to_netcdf(path = f'{mod[1]}/{filename}_{_date}.nc4', mode ='w',engine='scipy')
        var_final.close()


    # loadded = np.load(model_0 + f'/EDDI_{_date}.npy',allow_pickle=True)
    # #Add all 4 models to the same file
    # var_OUT.ETo[:,0,:,:,:] = np.load(model_0 + f'/EDDI_{_date}.npy',allow_pickle=True)
    # var_OUT.ETo[:,1,:,:,:] = np.load(model_1 + f'/EDDI_{_date}.npy',allow_pickle=True)
    # var_OUT.ETo[:,2,:,:,:] = np.load(model_2 + f'/EDDI_{_date}.npy',allow_pickle=True)
    # var_OUT.ETo[:,3,:,:,:] = np.load(model_3 + f'/EDDI_{_date}.npy',allow_pickle=True)
