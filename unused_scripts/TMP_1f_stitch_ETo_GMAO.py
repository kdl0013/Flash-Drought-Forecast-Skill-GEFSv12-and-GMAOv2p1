#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine together individual models from anomaly calculation and save into 
a single netcdf file for further processing.
"""

import xarray as xr
import numpy as np
import os
from glob import glob

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod = 'GMAO'
var = 'ETo'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#create dictionary of with 4 models (GMAO specifically)
model_dirs = {}

if var == 'EDDI':
    for mod in [0,1,2,3]:
        model_dirs[f'Model {mod}'] = f'{home_dir}/EDDI_mod{mod}/already_completed'
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


print(f'Saving stitched {var}.nc4 into {home_dir}.')
for _date in init_date_list:
    
    if var == 'EDDI':
        filename = 'EDDI'
        desc = 'Evaporative Demand Drought Index'
    elif var == 'RZSM':
        filename = f'{var}_anomaly'
        desc = 'RZSM anomaly. 7-day average mean lead time and 7-day window'
    elif var == 'ETo':
        filename = f'{var}_anomaly'
        desc = 'ETo anomaly. 7-day average mean by lead time and 7-day window.'
    
    #Base layer template to append to 
    template = xr.open_dataset(f'ETo_{_date}.nc4')
    var_OUT = xr.zeros_like(template)
        
    if var == 'EDDI':
        
        #Add data from each model into the empty file
        for mod in model_dirs.items():
            var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = np.load(mod[1] + f'/{filename}_{_date}.npy',allow_pickle=True)
            #print(mod)
            
        #save file as EDDI as netcdf
        var_final = xr.Dataset(
            data_vars = dict(
                EDDI = (['S', 'model','lead','Y','X'], var_OUT.ETo.values),
            ),
            coords = dict(
                X = var_OUT.X.values,
                Y = var_OUT.Y.values,
                lead = var_OUT.lead.values,
                model = var_OUT.model.values,
                S = template.S.values
            ),
            attrs = dict(
                Description = f'{desc}'),
        )                    
        
        #Save as a netcdf for later processing
        var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc4', mode ='w')
        
    elif var == 'RZSM':

        #Add data from each model into the empty file
        for mod in model_dirs.items():
            xr_file = xr.open_dataset(mod[1] + f'/{filename}_{_date}.nc4')
            
            try:
                var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = xr_file.Variable[:,:,:,:]
            except AttributeError:
                var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = xr_file.Variables[:,:,:,:]
            #print(mod)
        
        #ncview doesnt' like infinity values
        var_OUT_masked = var_OUT.where(abs(var_OUT.ETo) < 100)
        
        #save file as EDDI as netcdf
        var_final = xr.Dataset(
            data_vars = dict(
                RZSM_anom = (['S', 'model','lead','Y','X'], var_OUT_masked.ETo.values),
            ),
            coords = dict(
                X = var_OUT.X.values,
                Y = var_OUT.Y.values,
                lead = var_OUT.lead.values,
                model = var_OUT.model.values,
                S = template.S.values
            ),
            attrs = dict(
                Description = f'{desc}'),
        )                    
        
        #Save as a netcdf for later processing
        var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc4', mode ='w')
        
    elif var == 'ETo':
    
        #Add data from each model into the empty file
        for mod in model_dirs.items():
            xr_file = xr.open_dataset(mod[1] + f'/{filename}_{_date}.nc4')
            
            try:
                var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = xr_file.Variable[:,:,:,:]
            except AttributeError:
                var_OUT.ETo[:,int(mod[0][-1]),:,:,:] = xr_file.Variables[:,:,:,:]
            #print(mod)
        
        #ncview doesnt' like infinity values
        var_OUT_masked = var_OUT.where(abs(var_OUT.ETo) < 100)
        
        #save file as EDDI as netcdf
        var_final = xr.Dataset(
            data_vars = dict(
                ETo_anom = (['S', 'model','lead','Y','X'], var_OUT_masked.ETo.values),
            ),
            coords = dict(
                X = var_OUT.X.values,
                Y = var_OUT.Y.values,
                lead = var_OUT.lead.values,
                model = var_OUT.model.values,
                S = template.S.values
            ),
            attrs = dict(
                Description = f'{desc}'),
        )                    
        
        #Save as a netcdf for later processing
        var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc4', mode ='w')
           

