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
mod = 'model_name'
var = 'variables'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='GMAO'

home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'vas')  

#%%#SMERGE Soil Moisture
#SubX.dat depth file

#I brought this function out of loop, no need to repeat it more than once
T_FILE = xr.open_dataset(glob('ETo_1999-01-10.nc')[0])

def make_empty_nc_files(init_date_list,T_FILE, var):
    print(f'Making empty .nc files for {var}.')
    # count=0
    
    # for var in ['EDDI','RZSM','ETo']:
    for _date in init_date_list:
       
        model_dirs = {}
        # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
        if var == 'EDDI':
            filename = 'EDDI'
            desc = 'Evaporative Demand Drought Index'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                
        elif var == 'RZSM':
            filename = f'{var}_anomaly'
            desc = f'RZSM anomaly SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
        elif var == 'ETo':
            filename = f'{var}_anomaly'
            desc = f'ETo anomaly SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
        
        # try:
        #     np.load(f'{home_dir}/julian_lead_{_date}.npy')
        # except FileNotFoundError:
        #     np.save(f'{home_dir}/julian_lead_{_date}.npy',day_julian_b) 
            
            
        #don't remake empty files
        try:
            xr.open_dataset(f'{home_dir}/{filename}_{_date}.nc4')
            
        except FileNotFoundError:
            template = xr.open_dataset(f'ETo_{_date}.nc')
            #initialize empty file
            template = xr.zeros_like(template)

            if var == 'RZSM':
                var_final = xr.Dataset(
                    data_vars = dict(
                        RZSM_anom = (['S','model','lead','Y','X'], template.ETo.values),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = template.lead.values,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                
            elif var == 'EDDI':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        EDDI = (['S','model','lead','Y','X'], template.ETo.values),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = template.lead.values,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                

                
            elif var == 'ETo':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        ETo_anom = (['S','model','lead','Y','X'], template.ETo.values),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = template.lead.values,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                

            var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc')
            #compress so that when I re-write the file, it is quicker
            # '''But this doesn't work after the file is re-read and re-saved'''
            # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
            # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
            var_final.close()

   
#Create EDDI files           
make_empty_nc_files(init_date_list = init_date_list,T_FILE = T_FILE, var=var)
#%%'''Make empty files to keep track of anomaly calculations'''

#Make a completed list file for EDDI, add new names to next code block
new_eddi = f'EDDI_completed_nc_{mod}.txt'
new_eto_anom = f'ETo_completed_anomaly_nc_{mod}.txt'
new_rzsm = f'RZSM_completed_anomaly_nc_{mod}.txt'

#Create a new file for each index to keep track of what has been completed 
#since this is a pretty slow process
for name in [new_eddi,new_eto_anom,new_rzsm]:
    try:
        completed_dates = np.loadtxt(f'{script_dir}/{name}',dtype='str')         
    except OSError:
        os.system(f'touch {script_dir}/{name}')
        name_index = f'{name}'.split('_')[0]
        os.system(f'echo "Completed {name_index}" > {script_dir}/{name}')
