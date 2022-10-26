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
obs_dir = f'{dir1}/Data/MERRA2' #reference evapotranspiration
output_ETo_dir = f'{obs_dir}/ETo_SubX_values' #Refernce ET output directory
output_RZSM_dir = f'{obs_dir}/RZSM_SubX_values' 

os.system(f'mkdir -p {output_ETo_dir} {output_RZSM_dir}')

model_='NRL'
global evap
evap='Penman'

subX_dir = f'{dir1}/Data/SubX/{model_}'
n_processors = 10

os.chdir(subX_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob(f'ETo_{evap}*.nc4'))
date_list = [i[-14:-4] for i in date_list]


# _date = date_list[0]
# _date='1999-12-06'
#%% ETo gridMET
def anomaly_ETo_Priestley_gridMET_SubX_creation(_date) -> float:    
    file_name = f'ETo_SubX_anomaly_{evap}_{model_}_{_date}.nc4'
    
    # model_='EMC'
    # evap = 'Priestley'
    
    

    try:
        # xr.open_dataset(f'{output_ETo_dir}/{file_name}')
        xr.open_dataset('test.nc') #only to make new files
        print(f'Already completed date {_date}. Saved in {output_ETo_dir}')
    except FileNotFoundError:
        print(f'Working on initialized day {_date} to find MERRA values from SubX models, leads, & coordinates and saving data into {output_ETo_dir}.')

        #Open up anomaly file to get proper formatting
        sub_file = xr.open_dataset(f'{subX_dir}/ETo_{evap}_{model_}_{_date}.nc4')
        out_file = xr.zeros_like(sub_file)
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        obs_file = xr.open_dataset(f'{obs_dir}/ETo_anomaly_{evap}_MERRA.nc').astype('float64')
        '''Issue with the dates'''
        dates_list = list(obs_file.time.values)
        df = pd.DataFrame({'dates':dates_list})
        df['dates'] = pd.to_datetime(df['dates'])
        df['dates'] = df['dates'].dt.floor('d')
        
        new_date_list = np.array(df['dates'])
        
        obs_file=obs_file.assign_coords(time=new_date_list)
        
        #Mask for CONUS (don't process additional data)
        HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
        HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files
           
        '''Now that we have the subx reference ET file open, now we need to open up the 
        gridMET time series (all a single file). Select the same dates as SubX for 
        the full time series.
        
        Goal is to take each SubX file, and find the corresponding historical observed file (Y,X value), 
        and then append the SubX value and historical value to seperate arrays. Then 
        we can run the scatterplot to determine how close the observed is with SubX.
        
        #Next step is to make a copy of the subX file and just fill in that copied dataset
        #with values from gridMET gleam ref ET. This will assist with the ravel() function
        #to be used later
        
        Find the same date(s) values from each dataset and append to output dataset
        '''
     
        for i_lead in range(sub_file.ETo.shape[2]):
            #Because the S value is actually the first day of the forecast, don't add one
            date_val = pd.to_datetime(pd.to_datetime(_date) + dt.timedelta(days=i_lead))
            #1 day appears to have not been calculated because of julian days and 
            #anomaly code for +/-42 day (specically December 31, 2000)
            #Just take the average of the other two days before and after
            if np.count_nonzero( obs_file.ETo_anom.sel(time = date_val).values == 0) == 1593:
                date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                
                out_file.ETo[0,:, i_lead, :, :] = \
                np.nanmean( obs_file.ETo_anom.sel(time = date_val1).values + \
                     obs_file.ETo_anom.sel(time = date_val2).values)
            else:
                out_file.ETo[0,:, i_lead, :, :] = \
                    obs_file.ETo_anom.sel(time = date_val).values
                
                
        #Only keep 1 model
        out_file = out_file.isel(M=0)
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_anom = (['S','L','Y','X'], out_file.ETo.values),
            ),
            coords = dict(
                S = np.atleast_1d(pd.to_datetime(_date)),
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                L = np.arange(sub_file.L.shape[0]),
        
            ),
            attrs = dict(
                Description = f'MERRA2 reference ETo anomaly values on the exact same date and grid \
                cell as {model_} SubX data'),
        )                    
        

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')
        
                # '''Compress to save space since I do not have to write again to file'''
                # os.system(f'ncks -4 -L 1 {output_ETo_dir}/{file_name} {output_ETo_dir}/eto_test.nc4')
                # os.system(f'mv {output_ETo_dir}/eto_test.nc4 {output_ETo_dir}/{file_name}')

        print(f'Completed ETo_SubX_{_date}')  
        
        return()

#%%
if __name__ == '__main__':
    p = Pool(n_processors)
    p.map(anomaly_ETo_Priestley_gridMET_SubX_creation,date_list)
