#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Location of GMAO GEOS-5 hindcasts https://gmao.gsfc.nasa.gov/gmaoftp/gmaofcst/subx/GEOS_S2S_V2.1_fcst/IRI/

File name order: cape_GMAOGEOS_01apr2018_00z_d01_d45_m01.nc 
--- Only get every 5 days. Starts on 25jul2017.
---Need to get all 4 models

@author: kdl
"""

import xarray as xr
import numpy as np
import pandas as pd
import datetime as dt
import os
from glob import glob

dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
preprocess_dir = f'{dir1}/Data/SubX/GMAO/preprocess_2017_2021'
os.chdir(preprocess_dir)
out_dir = f'{dir1}/Data/SubX/GMAO'
vars_to_download = ['dswrf_sfc', 'uas_10m','vas_10m', 'tasmax_2m','tasmin_2m', 'mrso', 'spfh_2m','tdps_2m']

#
start_date = pd.to_datetime('07-29-2017')
end_date = pd.to_datetime('12-31-2021')

final_date_list = start_date + np.timedelta64(5,'D')

#dates
dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

dates_out = []
for i in dates:
    if i.month == 2 and i.day == 29:
        pass
    else:
        dates_out.append(i)

#Only every 5th date
dates=dates_out[::5]



#Create empty dataset to store outputs
var_OUT = np.zeros(shape=(1,4,45,27,59))


#name of file is weird, convert date to open file
def convert_date(_date):
    day = f'{_date.day:02}'
    year = _date.year
    month = _date.strftime('%B').lower()[0:3]
    out_date = f'{day}{month}{year}'
    return(out_date)
        
out_date = convert_date(start_date)
        
def output_file_name(var,out_date):
    varname = var.split('_')[0] + '_GMAO_'
    date_name = f'{pd.to_datetime(out_date).year}-{pd.to_datetime(out_date).month:02}-{pd.to_datetime(out_date).day:02}_m{model}.nc4'
    out_name = f'{varname}{date_name}'
    return(out_name)


#add to url with variable
for _date in dates:
    out_date = convert_date(_date)
    for var in vars_to_download:
        all_models = sorted(glob(f'{var}_GMAOGEOS_{out_date}*'))
        for mod_number,file in enumerate(all_models):
            
            file_o = xr.open_dataset(file)    
            var_name_file = list(file_o.data_vars)[0]
            var_OUT[:,mod_number,:,:,:] = file_o[var_name_file].values
            
        if 'tasmax' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    tasmax = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'Tmax GMAO GEOS-5. All ensembles and leads in one file'),
            )  
        elif 'tasmin' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    tasmin = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'Tmin GMAO GEOS-5. All ensembles and leads in one file'),
            )  
        elif 'vas' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    vas = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'V component wind GMAO GEOS-5. All ensembles and leads in one file'),
            )  
        elif 'uas' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    uas = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'U component wind GMAO GEOS-5. All ensembles and leads in one file'),
            )  
        elif 'mrso' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    mrso = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'Soil moisture GMAO GEOS-5. All ensembles and leads in one file'),
            )  
        elif 'tdps' in var:
            GMAO_out = xr.Dataset(
                data_vars = dict(
                    tdps = (['S','model','lead','Y','X'], var_OUT[:,:,:,:,:]),
                ),
                coords = dict(
                  
                    X = file_o.X.values,
                    Y = file_o.Y.values,
                    lead = range(45),
                    model = range(4),
        
                ),
                attrs = dict(
                    Description = 'Dewpoint temp GMAO GEOS-5. All ensembles and leads in one file'),
            )  
            
            out_name1 = 'out_1.nc4'
            GEFS_out.to_netcdf(path = f"{home_dir}/{var}/{out_name1}")
            #Now compress
            os.system(f'ncks -O -4 -L 1 {home_dir}/{var}/{out_name1} {save_out_dir_GEFS}/{final_out_name}')
            #Remove the old files
            os.system(f'rm {home_dir}/{var}/{out_name1}')
#%%

        
            out_name = output_file_name(var,out_date)
            # output.append(f'wget -nc -O {url1} {out_dir}/{out_name} &') #doesn't work
            output.append(f'wget -nc {url1} &')
            count+=1
                
    #Now only get data for every 5 days
    start_date = start_date + np.timedelta64(5,'D')
        
        

        
        
