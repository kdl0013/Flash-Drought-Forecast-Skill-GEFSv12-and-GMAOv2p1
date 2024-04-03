#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
I received 5 new dates for ESRL FIMr1p1. They are already pre-processed,
just change the names and move the files.

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob
import pandas as pd
import datetime as dt

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='RSMAS'
subX_dir = f'{dir1}/Data/SubX/ESRL/from_ESRL_agency_missing_days/keep_variables'
out_dir = f'{dir1}/Data/SubX/ESRL'

os.chdir(subX_dir)

#get file list
file_list_tas = sorted(glob('tas*.nc'))
file_list_rad = sorted(glob('rad*.nc'))

# zip_list = zip(file_list_tas,file_list_rad)

def change_date(file):
    split_ = file.split('_')
    date_ = split_[3]
    combine = f'{pd.to_datetime(date_).year}-{pd.to_datetime(date_).month:02}-{pd.to_datetime(date_).day:02}'
    
    name = f'{split_[0]}_ESRL_{combine}.nc4'
    return(name)

#TODO: Need to combine all 4 models into 1

date_unique = np.unique([i.split('_')[3] for i in file_list_tas])
#open file, resave elsewhere
for date_u in date_unique:
    tas_list = sorted(glob(f'ta*{date_u}*nc4'))
    template = xr.open_dataset(tas_list[0])
    template_o = np.empty((1,4,template.time.shape[0],template.Y.shape[0],template.X.shape[0]))
    for mod,file in enumerate(tas_list):
        file_o = xr.open_dataset(file)
        template_o[:,mod,:,:,:] = file_o[list(file_o.keys())[0]].values
    
    var_OUT = xr.Dataset(
        data_vars = dict(
            tas = (['S', 'M','L','Y','X'], template_o[:,:,:,:,:]),
        ),
        coords = dict(
            X = template.X.values,
            Y = template.Y.values,
            L = np.arange(len(template.time.values)),
            M = np.arange(4),
            
        ),
        attrs = dict(
            Description = 'Average temperature ESRL.',)
    )  

    out_name = change_date(file)
    var_OUT.to_netcdf(f'{out_dir}/{out_name}')
    
#open file, resave elsewhere
for date_u in date_unique:
    tas_list = sorted(glob(f'rad*{date_u}*nc4'))
    template = xr.open_dataset(tas_list[0])
    template_o = np.empty((1,4,template.time.shape[0],template.Y.shape[0],template.X.shape[0]))
    for mod,file in enumerate(tas_list):
        file_o = xr.open_dataset(file)
        template_o[:,mod,:,:,:] = file_o[list(file_o.keys())[0]].values
    
    var_OUT = xr.Dataset(
        data_vars = dict(
            rad = (['S', 'M','L','Y','X'], template_o[:,:,:,:,:]),
        ),
        coords = dict(
            X = template.X.values,
            Y = template.Y.values,
            L = np.arange(len(template.time.values)),
            M = np.arange(4),
            
        ),
        attrs = dict(
            Description = 'Net radiation.',)
        )
    out_name = change_date(file)
    var_OUT.to_netcdf(f'{out_dir}/{out_name}')
        

