#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will find create a new dataset of actual observations, but with
the format of SubX data. This will make correlation between datasets easier 
to calculate.
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from multiprocessing import Pool

#open a single file for SubX reference ET and a single file from gridMET 
dir1 = 'main_dir'
num_processors = int('procs')

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

subX_dir = f'{dir1}/Data/SubX/GMAO'
sub_split = subX_dir.split('/')[-1]
gridMET_dir = f'{dir1}/Data/gridMET' #reference evapotranspiration
output_ETo_dir = f'{gridMET_dir}/ETo_SubX_values' #Refernce ET output directory
smerge_in_dir = f'{dir1}/Data/SMERGE_SM/Raw_data' #raw files that have been boxed in CONUS
SM_SubX_out_dir = f'{dir1}/Data/SMERGE_SM/SM_SubX_values' #Smerge values overlayed on SubX grid
subx_anomaly_dir = f'{subX_dir}/anomaly'
#Didn't do EDDI because I may not need to do it.
# eddi_dir = f'{dir1}/Data/EDDI/convert_2_nc'
# eddi_subX_dir = f'{dir1}/Data/EDDI/EDDI_SubX_values'

image_dir = f'{dir1}/Outputs/MSE_plots'

#Make new directories
os.system(f'mkdir -p {image_dir}')
os.system(f'mkdir {SM_SubX_out_dir}')
# os.system(f'mkdir {eddi_subX_dir}')
os.system(f'mkdir {output_ETo_dir}')

os.chdir(subX_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob('mrso*_*-*.nc4'))
date_list = [i[-14:-4] for i in date_list]

# _date = date_list[0]
# _date='1999-12-06'
#%% ETo gridMET shortwave radiation (dswrf)
def gridMET_SubX_creation_dswrf(_date) -> float:    
    
    sub_name='dswrf'
    long_sub_name = 'surface_downwelling_shortwave_flux_in_air' #name in actual gridmet file
    gridName = 'srad_remap_final'
    try:
        xr.open_dataset(f'{output_ETo_dir}/{sub_name}_SubX_{_date}.nc4')
        print(f'Already completed date {_date} for {sub_name}.')
    except FileNotFoundError:
        
        
        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/{gridName}.nc').astype('float64')

        print(f'Working on initialized day {_date} for var {sub_name}.')
        sub_file = xr.open_dataset(f'{subX_dir}/dswrf_{sub_split}_{_date}.nc4')
        
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
        gridMET_out = (np.zeros_like(sub_file[f'{sub_name}'])).astype('float64')

        for i_lead in range(sub_file[f'{sub_name}'].shape[2]):
            #Need to add 1 to the final date_val value to be the correct date
            add_one = 1
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+add_one) 
     
            for i_Y in range(sub_file[f'{sub_name}'].shape[3]):
                for i_X in range(sub_file[f'{sub_name}'].shape[4]):
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero( gridMET_file[f'{long_sub_name}'].sel(day = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean( gridMET_file[f'{long_sub_name}'].sel(day = date_val1).isel(X = i_X, Y = i_Y).values + \
                                 gridMET_file[f'{long_sub_name}'].sel(day = date_val2).isel(X = i_X, Y = i_Y).values)
                        else:
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                                gridMET_file[f'{long_sub_name}'].sel(day = date_val).isel(X = i_X, Y = i_Y).values
                                     
        gridMET_get_weekly_leads = gridMET_out[:,:,::7,:,:]
        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                dswrf = (['S','model','lead','Y','X'], gridMET_get_weekly_leads[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = np.arange(0,7),
                model = sub_file.M.values,
        
            ),
            attrs = dict(
                Description = f'gridMET shortwave radiation {sub_name} values on the exact same date and grid \
                cell as SubX data'),
        )                    
        
        file_name = f'{sub_name}_SubX_{_date}.nc4'

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')

        
        # '''Compress to save space since I do not have to write again to file'''
        # os.system(f'ncks -4 -L 1 {output_ETo_dir}/{file_name} {output_ETo_dir}/eto_test.nc4')
        # os.system(f'mv {output_ETo_dir}/eto_test.nc4 {output_ETo_dir}/{file_name}')

        print(f'Completed {sub_name}_SubX_{_date}.')
#%% (tasmax)
def gridMET_SubX_creation_tasmax(_date) -> float:    
    
    sub_name='tasmax'
    gridName = 'tmmx_remap_final'
    long_sub_name = 'air_temperature' #name in actual gridmet file

    try:
        xr.open_dataset(f'{output_ETo_dir}/{sub_name}_SubX_{_date}.nc4')
        print(f'Already completed date {_date} for {sub_name}.')
    except FileNotFoundError:
        
        
        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/{gridName}.nc').astype('float64')

        print(f'Working on initialized day {_date} for var {sub_name}.')
        sub_file = xr.open_dataset(f'{subX_dir}/{sub_name}_{sub_split}_{_date}.nc4')
        
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
        gridMET_out = (np.zeros_like(sub_file[f'{sub_name}'])).astype('float64')

        for i_lead in range(sub_file[f'{sub_name}'].shape[2]):
            #Need to add 1 to the final date_val value to be the correct date
            add_one = 1
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+add_one) 
     
            for i_Y in range(sub_file[f'{sub_name}'].shape[3]):
                for i_X in range(sub_file[f'{sub_name}'].shape[4]):
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero( gridMET_file[f'{long_sub_name}'].sel(day = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean( gridMET_file[f'{long_sub_name}'].sel(day = date_val1).isel(X = i_X, Y = i_Y).values + \
                                 gridMET_file[f'{long_sub_name}'].sel(day = date_val2).isel(X = i_X, Y = i_Y).values)
                        else:
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                                gridMET_file[f'{long_sub_name}'].sel(day = date_val).isel(X = i_X, Y = i_Y).values
                                     
        gridMET_get_weekly_leads = gridMET_out[:,:,::7,:,:]
        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                tasmax = (['S','model','lead','Y','X'], gridMET_get_weekly_leads[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = np.arange(0,7),
                model = sub_file.M.values,
        
            ),
            attrs = dict(
                Description = f'gridMET max temp {sub_name} values on the exact same date and grid \
                cell as SubX data'),
        )                    
        
        file_name = f'{sub_name}_SubX_{_date}.nc4'

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')

#%% (tasmin)
def gridMET_SubX_creation_tasmin(_date) -> float:    
    
    sub_name='tasmin'
    gridName = 'tmmn_remap_final'
    long_sub_name = 'air_temperature' #name in actual gridmet file

    try:
        xr.open_dataset(f'{output_ETo_dir}/{sub_name}_SubX_{_date}.nc4')
        print(f'Already completed date {_date} for {sub_name}')
    except FileNotFoundError:
        
        
        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/{gridName}.nc').astype('float64')

        print(f'Working on initialized day {_date} for var {sub_name}.')
        sub_file = xr.open_dataset(f'{subX_dir}/{sub_name}_{sub_split}_{_date}.nc4')
        
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
        gridMET_out = (np.zeros_like(sub_file[f'{sub_name}'])).astype('float64')

        for i_lead in range(sub_file[f'{sub_name}'].shape[2]):
            #Need to add 1 to the final date_val value to be the correct date
            add_one = 1
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+add_one) 
     
            for i_Y in range(sub_file[f'{sub_name}'].shape[3]):
                for i_X in range(sub_file[f'{sub_name}'].shape[4]):
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero( gridMET_file[f'{long_sub_name}'].sel(day = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean( gridMET_file[f'{long_sub_name}'].sel(day = date_val1).isel(X = i_X, Y = i_Y).values + \
                                 gridMET_file[f'{long_sub_name}'].sel(day = date_val2).isel(X = i_X, Y = i_Y).values)
                        else:
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                                gridMET_file[f'{long_sub_name}'].sel(day = date_val).isel(X = i_X, Y = i_Y).values
                                     
        gridMET_get_weekly_leads = gridMET_out[:,:,::7,:,:]
        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                tasmin = (['S','model','lead','Y','X'], gridMET_get_weekly_leads[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = np.arange(0,7),
                model = sub_file.M.values,
        
            ),
            attrs = dict(
                Description = f'gridMET min temp {sub_name} values on the exact same date and grid \
                cell as SubX data'),
        )                    
        
        file_name = f'{sub_name}_SubX_{_date}.nc4'

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')
        

#%%
def gridMET_SubX_creation_humidity(_date) -> float:    
    
# for _date in date_list:
    sub_name='RelativeHumidity'
    out_name = 'RH'
    gridName = 'RelativeHumidity_gridMET_merged'
    long_sub_name = 'RH' #name in actual gridmet file

    try:
        xr.open_dataset(f'{output_ETo_dir}/{out_name}_SubX_{_date}.nc4')
        print(f'Already completed date {_date} for {sub_name}')
    except FileNotFoundError:

        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/{gridName}.nc').astype('float64')

        print(f'Working on initialized day {_date} for var {sub_name}.')
        sub_file = xr.open_dataset(f'{subX_dir}/{sub_name}_{_date}.nc4')
  
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
        gridMET_out = (np.zeros_like(sub_file[f'{long_sub_name}'])).astype('float64')

        for i_lead in range(sub_file[f'{long_sub_name}'].shape[2]):
            #Need to add 1 to the final date_val value to be the correct date
            add_one = 1
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+add_one) 
        
            for i_Y in range(sub_file[f'{long_sub_name}'].shape[3]):
                for i_X in range(sub_file[f'{long_sub_name}'].shape[4]):
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero( gridMET_file[f'{long_sub_name}'].sel(day = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean( gridMET_file[f'{long_sub_name}'].sel(day = date_val1).isel(lon = i_X, lat = i_Y).values + \
                                 gridMET_file[f'{long_sub_name}'].sel(day = date_val2).isel(lon = i_X, lat = i_Y).values)
                        else:
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                                gridMET_file[f'{long_sub_name}'].sel(day = date_val).isel(lon = i_X, lat = i_Y).values
                                     
        gridMET_get_weekly_leads = gridMET_out[:,:,::7,:,:]
        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                RH = (['S','model','lead','Y','X'], gridMET_get_weekly_leads[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = np.arange(0,7),
                model = sub_file.model.values,
        
            ),
            attrs = dict(
                Description = f'gridMET relative humidity {long_sub_name} values on the exact same date and grid \
                cell as SubX data'),
        )                    
        
        file_name = f'{out_name}_SubX_{_date}.nc4'

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')
#%%Run all four functions
if __name__ == '__main__':
    p = Pool(num_processors)
    print('Working on converting shortwave radiation (dswrf) to SubX format to assist with anomaly correlation')
    p.map(gridMET_SubX_creation_dswrf, date_list)
    print('Working on maximum temperature (tasmax) to SubX format to assist with anomaly correlation')
    p.map(gridMET_SubX_creation_tasmax, date_list)
    print('Working on minimim temperature (tasmin) to SubX format to assist with anomaly correlation')
    p.map(gridMET_SubX_creation_tasmin, date_list)
    print('Working on relative humidity to SubX format to assist with anomaly correlation')
    p.map(gridMET_SubX_creation_humidity, date_list)
