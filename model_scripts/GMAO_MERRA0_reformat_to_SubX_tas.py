#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will find create a new dataset of actual observations, but with
the format of SubX data. This will make correlation between datasets easier 
to calculate.

FOR GMAO GEOS-5, file name is technically day 0.

"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from multiprocessing import Pool




dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod='GMAO'
var = 'tas'
n_processors = 8
home_dir = f'{dir1}/Data/SubX/{mod}'
obs_dir = f'{dir1}/Data/MERRA2' #reference evapotranspiration

#%%

#CONUS mask
conus_mask = xr.open_dataset(f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4')

os.chdir(home_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob(f'tasmax_{mod}*.nc4'))
date_list = [i[-14:-4] for i in date_list]


# _date = date_list[0]
# _date='1999-12-06'
#%%
def convert_OBS_to_SubX(_date) -> float:   
    # var='tas'
    # date_list = sorted(glob(f'{var}_{mod}*.nc4'))
    # date_list = [i[-14:-4] for i in date_list]
    if var == 'srad':
        merra = xr.open_dataset(f'{obs_dir}/radiation_anomaly_MERRA.nc')
    else:
        merra = xr.open_dataset(f'{obs_dir}/{var}_anomaly_MERRA.nc')

    output_obs_dir = f'{obs_dir}/{var}_SubX_values' #Refernce ET output directory
    os.system(f'mkdir -p {output_obs_dir}')
    obs_file_name = f'{var}_SubX_{mod}_{_date}.nc4'
    # os.system(f'mkdir -p {output_obs_dir}')
    
    try:
        xr.open_dataset(f'{output_obs_dir}/{obs_file_name}')
        # xr.open_dataset('test.nc') #only to make new files
        print(f'Already completed date {_date}. Saved in {output_obs_dir}')
    except FileNotFoundError:
        print(f'Working on initialized day {_date} to find MERRA values from SubX models, leads, & coordinates and saving data into {output_obs_dir}.')
        
        
        #Open up anomaly file to get proper formatting
        sub_file = xr.open_dataset(f'{home_dir}/{var}_{mod}_{_date}.nc4')
        out_file = xr.zeros_like(sub_file)
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        '''Issue with the dates'''
        dates_list = list(merra.time.values)
        df = pd.DataFrame({'dates':dates_list})
        df['dates'] = pd.to_datetime(df['dates'])
        df['dates'] = df['dates'].dt.floor('d')
        
        new_date_list = np.array(df['dates'])
        
        merra=merra.assign_coords(time=new_date_list)
        
        output_f = xr.zeros_like(sub_file)
 
        for i_lead in range(sub_file[list(sub_file.keys())[0]].L.shape[0]):
            #Because the S value is actually the first day of the forecast, don't add one
            date_val = pd.to_datetime(pd.to_datetime(_date) + dt.timedelta(days=i_lead))
            #1 day appears to have not been calculated because of julian days and 
            #anomaly code for +/-42 day (specically December 31, 2000)
            #Just take the average of the other two days before and after
            if np.count_nonzero(merra[list(merra.keys())[0]].sel(time = date_val).values == 0) == 1593:
                date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                
                output_f[list(output_f.keys())[0]][0,:, i_lead, :, :] = \
                np.nanmean( merra[list(merra.keys())[0]].sel(time = date_val1).values + \
                     merra[list(merra.keys())[0]].sel(time = date_val2).values)
            else:
                output_f[list(output_f.keys())[0]][0,:, i_lead, :, :] = \
                    merra[list(merra.keys())[0]].sel(time = date_val).values
                
                
        #Only keep 1 model
        out_file = output_f.isel(M=0)
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                temp_mean = (['S','L','Y','X'],   out_file[list(output_f.keys())[0]].values),
            ),
            coords = dict(
                S = np.atleast_1d(pd.to_datetime(_date)),
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                L = np.arange(sub_file.L.shape[0]),
        
            ),
            attrs = dict(
                Description = f'MERRA2 {var} values on the exact same date and grid \
                cell as {mod} SubX data'),
        )                    
        
 
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_obs_dir}/{obs_file_name}', mode ='w')
        
                # '''Compress to save space since I do not have to write again to file'''
                # os.system(f'ncks -4 -L 1 {output_ETo_dir}/{file_name} {output_ETo_dir}/eto_test.nc4')
                # os.system(f'mv {output_ETo_dir}/eto_test.nc4 {output_ETo_dir}/{file_name}')
 
        # print(f'Completed {var}_SubX_{_date}')  
        
     

    return(0)
#%%
if __name__ == '__main__':
    p = Pool(n_processors)
    p.map(convert_OBS_to_SubX,date_list)
#%%
