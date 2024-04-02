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

gridMET_dir = f'{dir1}/Data/gridMET' #reference evapotranspiration
output_ETo_dir = f'{gridMET_dir}/ETo_SubX_values' #Refernce ET output directory
smerge_in_dir = f'{dir1}/Data/GLEAM_RZSM' #raw files that have been boxed in CONUS
SM_SubX_out_dir = f'{smerge_in_dir}/SM_SubX_values' #Smerge values overlayed on SubX grid
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
#%% ETo gridMET
def anomaly_ETo_gridMET_SubX_creation(_date) -> float:    
    try:
        xr.open_dataset(f'{output_ETo_dir}/ETo_SubX_anomaly_{_date}.nc4')
        print(f'Already completed date {_date}. Saved in {output_ETo_dir}')
    except FileNotFoundError:
        #Open up anomaly file to get proper formatting
        sub_anomaly = xr.open_dataset(f'{subx_anomaly_dir}/ETo_anomaly_{_date}.nc4')
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/ETo_anomaly_gridMET_merged.nc').astype('float64')

        print(f'Working on initialized day {_date} to find gridMET values from SubX models, leads, & coordinates and saving data into {output_ETo_dir}.')
        sub_file = xr.open_dataset(f'ETo_{_date}.nc')
        
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
        gridMET_out = (np.zeros_like(sub_file.ETo)).astype('float64')
     
        for i_lead in range(sub_file.ETo.shape[2]):
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
        
            for i_Y in range(sub_file.ETo.shape[3]):
                for i_X in range(sub_file.ETo.shape[4]):
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero( gridMET_file.ETo_gridmet.sel(day = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean( gridMET_file.ETo_gridmet.sel(day = date_val1).isel(lon = i_X, lat = i_Y).values + \
                                 gridMET_file.ETo_gridmet.sel(day = date_val2).isel(lon = i_X, lat = i_Y).values)
                        else:
                            gridMET_out[0,:, i_lead, i_Y, i_X] = \
                                gridMET_file.ETo_gridmet.sel(day = date_val).isel(lon = i_X, lat = i_Y).values
                                     
        gridMET_get_weekly_leads = gridMET_out[:,:,::7,:,:]
        #Create a new dataset in which we can append the average leads (3-4,3-5,3-6,4-6)
        #for anomalies
        final_lead_file = np.empty(shape=[gridMET_out.shape[0],gridMET_out.shape[1],sub_anomaly.ETo_anom.shape[2],gridMET_out.shape[3],gridMET_out.shape[4]])
        final_lead_file[:,:,0:7,:,:] = gridMET_get_weekly_leads[:,:,:,:,:]
        
        for Y in range(final_lead_file.shape[3]):
            for X in range(final_lead_file.shape[4]):
                if (HP_conus_mask.High_Plains[0,Y,X].values in np.arange(1,7)):
                    for idx,lead in enumerate(sub_anomaly.lead.values):
                        #Don't have to work on the first 7 values (they are already in the file)
                        if lead ==0 or lead ==1 or lead==2 or lead==3 or lead==4 or lead ==5 or lead==6:
                            pass
                        else:
                            if lead == 3.4:
                                lead1 = 3
                                lead2 = 4
                                #Take averages, add to file
                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X]])
                            
                            elif lead == 3.5:
                                
                                lead1 = 3
                                lead2 = 4
                                lead3 = 5

                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X]])

                            elif lead == 3.6:
                                
                                lead1 = 3
                                lead2 = 4
                                lead3 = 5
                                lead4 = 6

                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X],final_lead_file[0,:,lead4,Y,X]])

                            elif lead == 4.6:
                                
                                lead1 = 4
                                lead2 = 5
                                lead3 = 6

                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X]])

        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_anom = (['S','model','lead','Y','X'], final_lead_file[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = sub_anomaly.lead.values,
                model = sub_file.model.values,
        
            ),
            attrs = dict(
                Description = 'gridMET reference ETo anomaly values on the exact same date and grid \
                cell as SubX data'),
        )                    
        
        file_name = f'ETo_SubX_anomaly_{_date}.nc4'

        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/{file_name}', mode ='w')
        
        # '''Compress to save space since I do not have to write again to file'''
        # os.system(f'ncks -4 -L 1 {output_ETo_dir}/{file_name} {output_ETo_dir}/eto_test.nc4')
        # os.system(f'mv {output_ETo_dir}/eto_test.nc4 {output_ETo_dir}/{file_name}')

        print(f'Completed ETo_SubX_{_date}')


#%%   
def anomaly_GLEAM_SubX_creation(_date):    
    try:
        xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_anomaly_{_date}.nc4')
        print(f'Already completed date {_date}. Saved in {SM_SubX_out_dir}.')
        
    except FileNotFoundError:
        #Open up anomaly file to get proper formatting
        sub_anomaly = xr.open_dataset(f'{subx_anomaly_dir}/RZSM_anomaly_{_date}.nc4')
        # _date = ''
        print(f'Working on initialized day {_date} to find GLEAM anomaly values from SubX models, leads, & coordinates and saving data into {SM_SubX_out_dir}.')
        sub_file = xr.open_dataset(f'mrso_GMAO_{_date}.nc4')
        
        SMERGE_file = xr.open_dataset(f'{smerge_in_dir}/RZSM_anomaly_GLEAM.nc')
        
        '''Same format as SubX. Find the correct dates, leads, models and just fill in'''
        smerge_out = (np.zeros_like(sub_file.mrso)).astype('float64')
        
        #Mask for CONUS (don't process additional data)
        HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
        HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files
        
        #add the values from gleam 
        for i_lead in range(sub_file.mrso.shape[2]):
            date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
        
            for i_Y in range(sub_file.mrso.shape[3]):
                for i_X in range(sub_file.mrso.shape[4]):
                    #high plains mask is only for CONUS
                    if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                    
                        #1 day appears to have not been calculated because of julian days and 
                        #anomaly code for +/-42 day (specically December 31, 2000)
                        #Just take the average of the other two days before and after
                        if np.count_nonzero(SMERGE_file.RZSM.sel(time = date_val).values == 0) == 1593:
                            date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                            date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                            
                            smerge_out[0,:, i_lead, i_Y, i_X] = \
                            np.nanmean(SMERGE_file.RZSM.sel(time = date_val1).isel(lon = i_X, lat = i_Y).values + \
                                SMERGE_file.RZSM.sel(time = date_val2).isel(lon = i_X, lat = i_Y).values)
                        else:
                            smerge_out[0,:, i_lead, i_Y, i_X] = \
                                SMERGE_file.RZSM.sel(time = date_val).isel(lon = i_X, lat = i_Y).values
        
    
        smerge_get_weekly_leads = smerge_out[:,:,::7,:,:]
        #Create a new dataset in which we can append the average leads (3-4,3-5,3-6,4-6)
        #for anomalies
        final_lead_file = np.empty(shape=[smerge_out.shape[0],smerge_out.shape[1],sub_anomaly.RZSM_anom.shape[2],smerge_out.shape[3],smerge_out.shape[4]])
        final_lead_file[:,:,0:7,:,:] = smerge_get_weekly_leads[:,:,:,:,:]
        
        #Create weekly leads all in one file
        '''We have created weekly leads for several weeks (see below)
        out_lead = np.array([0,1,2,3,4,5,6,3.4,3.5,3.6,4.6])
        
        Weeks: 0,1,2,3,4. Then we take averages of multiple weeks including:
            weeks 3&4, 3&4&5, 3&4&5&6, 4&5&6.

        '''

        for Y in range(final_lead_file.shape[3]):
            for X in range(final_lead_file.shape[4]):
                if (HP_conus_mask.High_Plains[0,Y,X].values in np.arange(1,7)):
                    for idx,lead in enumerate(sub_anomaly.lead.values):
                        #Don't have to work on the first 7 values (they are already in the file)
                        if lead ==0 or lead ==1 or lead==2 or lead==3 or lead==4 or lead ==5 or lead==6:
                            pass
                        else:
                            if lead == 3.4:
                                lead1 = 3
                                lead2 = 4
                                
                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X]])
                            
                            elif lead == 3.5:
                                
                                lead1 = 3
                                lead2 = 4
                                lead3 = 5

                                final_lead_file[0,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X]])

                            elif lead == 3.6:
                                
                                lead1 = 3
                                lead2 = 4
                                lead3 = 5
                                lead4 = 6

                                final_lead_file[0,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X],final_lead_file[0,:,lead4,Y,X]])

                            elif lead == 4.6:
                                
                                lead1 = 4
                                lead2 = 5
                                lead3 = 6

                                final_lead_file[:,:,idx,Y,X] = \
                                    np.nanmean([final_lead_file[0,:,lead1,Y,X],final_lead_file[0,:,lead2,Y,X],\
                                               final_lead_file[0,:,lead3,Y,X]])

        
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                RZSM_anom = (['S','model','lead','Y','X'], final_lead_file[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = sub_anomaly.lead.values,
                model = sub_file.M.values,
        
            ),
            attrs = dict(
                Description = 'GLEAM RZSM (0-100cm) anomaly values (from m3/m3) on the exact same date and grid \
                cell as SubX data'),
        ) 

        #Save as a netcdf for later processing
        file_name = f'SM_SubX_anomaly_{_date}.nc4'
        var_OUT.to_netcdf(path = f'{SM_SubX_out_dir}/{file_name}', mode ='w')
        
        # '''Compress to save space since I do not have to write again to file'''
        # os.system(f'ncks -4 -L 1 {SM_SubX_out_dir}/{file_name} {SM_SubX_out_dir}/f_test.nc4')
        # os.system(f'mv {SM_SubX_out_dir}/f_test.nc4 {SM_SubX_out_dir}/{file_name}')

        #To find the correct date
        print(f'Completed SM_SubX_anomaly_{_date}')        


#%% EDDI (don't run...yet) 
# def EDDI_SubX_creation(_date):    
#     try:
#         xr.open_dataset(f'{eddi_subX_dir}/EDDI_SubX_{_date}.nc')
#         print(f'Already completed date {_date}. Saved in {eddi_subX_dir}')
#     except FileNotFoundError:
        
    
#         #Open gridMET and subset to the correct dates as SubX for full time series
#         eddi_file = xr.open_dataset(f'{eddi_dir}/eddi_merged.nc4').astype('float64')
#         eddi_file = eddi_file.sel(time=slice("1999-01-10","2016-02-09"))

#         print(f'Working on initialized day {_date} to retrieve EDDI values from SubX models, leads, & coordinates and saving data into {eddi_subX_dir}.')
#         var_file = xr.open_dataset(f'{subX_dir}/EDDI_{_date}.nc4')
        
#         #Remove dimension S[1] because it has not data (it's a GMAO GEOS specific issue)
#         var_file_arr = var_file.EDDI[:,:,:,:,:].values
        
        
#         '''Now that we have the subx reference ET file open, now we need to open up the 
#         gridMET time series (all a single file). Select the same dates as SubX for 
#         the full time series.
        
#         Goal is to take each SubX file, and find the corresponding historical observed file (Y,X value), 
#         and then append the SubX value and historical value to seperate arrays. Then 
#         we can run the scatterplot to determine how close the observed is with SubX.
        
#         #Next step is to make a copy of the subX file and just fill in that copied dataset
#         #with values from gridMET gleam ref ET. This will assist with the ravel() function
#         #to be used later
        
#         Find the same date(s) values from each dataset and append to output dataset
#         '''
#         eddi_out = (np.zeros_like(var_file_arr).squeeze()).astype('float64')
     
#         for model in range(var_file.EDDI.shape[1]):
#             for i_lead in range(var_file.EDDI.shape[2]):
#                 date_val = pd.to_datetime(var_file.S.values[0]) + dt.timedelta(days=i_lead)
            
#                 for i_Y in range(var_file.EDDI.shape[3]):
#                     for i_X in range(var_file.EDDI.shape[4]):
                        
#                         eddi_out[0,model, i_lead, i_Y, i_X] = \
#                             eddi_file.EDDI.sel(time = date_val).isel(X = i_X, Y = i_Y).values
        
#         #Convert to an xarray object
#         var_OUT = xr.Dataset(
#             data_vars = dict(
#                 ETo_SubX_value = (['S','model','lead','Y','X'], eddi_out[:,:,:,:,:]),
#             ),
#             coords = dict(
#                 S = var_file.S.values,
#                 X = var_file.X.values,
#                 Y = var_file.Y.values,
#                 lead = var_file.lead.values,
#                 model = var_file.model.values,
        
#             ),
#             attrs = dict(
#                 Description = 'Evaporative Demand Drought Index values on the exact same date and grid \
#                 cell as SubX data'),
#         )                    
            
#         #Save as a netcdf for later processing
#         var_OUT.to_netcdf(path = f'{eddi_subX_dir}/ETo_SubX_{_date}.nc4', mode ='w')
#         #To find the correct date
        
#         print(f'Completed ETo_SubX_{_date}')
# print(f'Working from directory {subX_dir}')


#%%Run all four functions
if __name__ == '__main__':
    p = Pool(num_processors)
    print('Working on converting SMERGE anomalies to SubX format to assist with anomaly correlation')
    p.map(anomaly_GLEAM_SubX_creation, date_list)
    print('Completed.')
    print('Working on converting ETo anomalies to SubX format to assist with anomaly correlation')    
    p.map(anomaly_ETo_gridMET_SubX_creation, date_list)
    # print('Working on EDDI.')
    # p.map(EDDI_SubX_creation, date_list)

#%%


# #%% #I don't need this code anymore, this isn't the correct correlation methodology
# '''Scatterplots for SMERGE RZSM and gridMET Reference ET compared to SubX 
# corresponding days and coordinates'''

# #Need to open each day, np.ravel the file, append to an empty array.
# #Once all files are opened, then produce a scatter plot

# def flatten_RZSM_ETo_EDDI(date_list):
#     '''Open SMERGE and SubX file'''
    
#     ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().values.ravel())    
#     file_length = len(date_list)
#     final_length = ravel_length * file_length
    
    
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_RZSM = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open SMERGE file and ravel, and append to empty array
#         smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].values.ravel()
#         #append smerge to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_RZSM[count_value:ravel_end,0] = smerge_file[:]
#         correlation_RZSM[count_value:ravel_end,1] = sub_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
        
#     '''Open gridMET and SubX file'''
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_ETo = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         g_file = xr.open_dataset(f'{output_dir}/ETo_SubX_{date_}.nc').to_array().values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         s_file = xr.open_dataset(f'{subX_dir}/ETo_{date_}.nc4').ETo[0,:,:,:,:].values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_ETo[count_value:ravel_end,0] = g_file[:]
#         correlation_ETo[count_value:ravel_end,1] = s_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
        
#     '''Open EDDI and SubX file'''
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_EDDI = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         g_file = xr.open_dataset(f'{output_dir}/EDDI_SubX_{date_}.nc').to_array().values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         s_file = xr.open_dataset(f'{subX_dir}/EDDI_{date_}.nc4').ETo[0,:,:,:,:].values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_EDDI[count_value:ravel_end,0] = g_file[:]
#         correlation_EDDI[count_value:ravel_end,1] = s_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
           
#     return(correlation_RZSM, correlation_ETo, correlation_EDDI)


# #TODO: keep working heres
# def average_realization(date_list):
#     '''Open gridMET and SubX file'''
#     ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().mean(dim='model').values.ravel())    
#     file_length = len(date_list)
    
#     final_length = ravel_length * file_length
    
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_RZSM_realization_mean = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         #Open SMERGE file and ravel, and append to empty array
#         smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().mean(dim='model').values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].mean(dim='model').values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_RZSM_realization_mean[count_value:ravel_end,0] = smerge_file[:]
#         correlation_RZSM_realization_mean[count_value:ravel_end,1] = sub_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
        
#     '''Open gridMET and SubX file'''
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_ETo_realization_mean = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         g_file = xr.open_dataset(f'{output_dir}/ETo_SubX_{date_}.nc').to_array().mean(dim='model').values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         s_file = xr.open_dataset(f'{subX_dir}/ETo_{date_}.nc4').ETo[0,:,:,:,:].mean(dim='model').values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_ETo_realization_mean[count_value:ravel_end,0] = g_file[:]
#         correlation_ETo_realization_mean[count_value:ravel_end,1] = s_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
        
#     '''Open EDDI and SubX file'''
#     ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().mean(dim='model').values.ravel())    
#     file_length = len(date_list)
    
#     final_length = ravel_length * file_length
    
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_EDDI_realization_mean = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         g_file = xr.open_dataset(f'{output_dir}/EDDI_SubX_{date_}.nc').to_array().mean(dim='model').values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         s_file = xr.open_dataset(f'{subX_dir}/EDDI_{date_}.nc4').ETo[0,:,:,:,:].mean(dim='model').values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_EDDI_realization_mean[count_value:ravel_end,0] = g_file[:]
#         correlation_EDDI_realization_mean[count_value:ravel_end,1] = s_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
#     return(correlation_RZSM_realization_mean,correlation_ETo_realization_mean,correlation_EDDI_realization_mean)





# '''Find only the realization mean and plot scatterplot and heatmap'''




# def plot_setup(date_list):
#     correlation_RZSM, correlation_ETo, correlation_EDDI = flatten_RZSM_ETo_EDDI(date_list)    
#     correlation_RZSM_realization_mean,correlation_ETo_realization_mean,correlation_EDDI_realization_mean = \
#         (average_realization(date_list))
    
#     #Find the areas where there are only values in the gridMET dataset because
#     #gridMET is already restricted to CONUS.
#     #make inf values np.nan
#     correlation_ETo[correlation_ETo > 1e308] = np.nan
#     #Make all negative values equal to 0 because Ref. ET cannot be negative
#     correlation_ETo[correlation_ETo < 0] = 0
    
#     no_nan = correlation_ETo[(~np.isnan(correlation_ETo[:,0]) & (~np.isnan(correlation_ETo[:,1]))) ]
    
#     #Same with mean
#     correlation_ETo_realization_mean[correlation_ETo_realization_mean > 1e308] = np.nan
#     #Make all negative values equal to 0 because Ref. ET cannot be negative
#     correlation_ETo_realization_mean[correlation_ETo_realization_mean < 0] = 0
    
#     no_nan_mean = correlation_ETo_realization_mean[(~np.isnan(correlation_ETo_realization_mean[:,0]) & (~np.isnan(correlation_ETo_realization_mean[:,1]))) ]
        
    
#     # #remove gridMET grid cells with np.nan because of outside CONUS bbox
#     # corr_1 = correlation_ETo[~np.isnan(correlation_ETo[:,0])]
#     # #remove SubX values with np.nan because of differences in grid coordinates/missing values in a few CONUS grid cells
#     # corr_2 = corr_1[~np.isnan(corr_1[:,1])]
    
#     # # Scatterplot and Correlations
#     # # Data
    
#     #GridMET is x
#     x = no_nan[:,0]
#     #SubX is y
#     y = no_nan[:,1]
    
#     #Mean realization
#     #GridMET is x    return(correlation_ETo_realization_mean)
#     mean_x = no_nan_mean[:,0]
#     #SubX is y
#     mean_y = no_nan_mean[:,1]        
    
#     return(x,y,mean_x,mean_y)


# #%%
# var = 'ETo'
# comparison_var = 'gridMET'
# comparison_type = 'mm/d'
# max_val = 30

# def plot_heatmaps_MSE(var,comparison_var,comparison_type,max_axis_val):
#     try:
#         plt.imread(f'{image_dir}/{var}_MSE.tif')
#         plt.imread(f'{image_dir}/{var}_heatmap_1.tif')
#         plt.imread(f'{image_dir}/{var}_heatmap_2.tif')
#         plt.imread(f'{image_dir}/{var}_MSE_mean.tif')
#     except FileNotFoundError:
#         print(f'Finding MSE and heatmap for {var}.')    
#         x,y,mean_x,mean_y = ETo_plot_setup(date_list)
        
#         x_line = np.linspace(0,max_axis_val,100)
#         y_line = np.linspace(0,max_axis_val,100)
        
        
        

    

#     try:
#         plt.imread(f'{image_dir}/RZSM_MSE_mean.tif')
#     except FileNotFoundError:
#         ax = sns.scatterplot(x = mean_a, y = mean_b,s=0.25)
#         ax.set_title("SubX RZSM vs. SMERGE RZSM")
#         ax.set_xlabel("SMERGE (m3/m3)")
#         ax.set_ylabel("SubX (m3/m3)")
#         ax.set(xlim=(0,max_RZSM_val))
#         ax.set(ylim=(0,max_RZSM_val))
#         plt.plot(x_line,y_line,color='r')
#         ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(mean_a,mean_b),5)}", fontsize=9) #add text
#         plt.savefig(f'{image_dir}/RZSM_MSE_mean.tif', dpi=300)
#         plt.close()
        
         
#     try:
#         plt.imread(f'{image_dir}/RZSM_heatmap_mean.tif')
#     except FileNotFoundError:
#     #Plot heatmap
#         heatmap, xedges, yedges = np.histogram2d(mean_a, mean_b, bins=100)
#         extent = [0, yedges[-1], 0, yedges[-1]]
#         plt.imshow(heatmap.T, extent=extent, origin='lower')
#         plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(mean_a,mean_b),5)}", fontsize=14, color = 'red')
#         plt.xlabel('SMERGE RZSM (m3/m3)')
#         plt.ylabel('SubX RZSM (m3/m3)')
#         plt.savefig(f'{image_dir}/RZSM_heatmap_mean.tif',dpi=300)
#         plt.close()        
        
#         if var =='ETO':
#             v1 = 1
#             v2 = 27
#             v3 = 3
#             v4= 20
#             v5= 12
#         elif var == 'RZSM':
#             v1 = 0.01
#             v2 = 0.94
#             v3 = 0.10
#             v4= 0.75
#             v5= 0.75

        
        
        
#         try:
#             plt.imread(f'{image_dir}/{var}_MSE.tif')
#         except FileNotFoundError:
#             ax = sns.scatterplot(x = x, y = y,s=0.25)
#             ax.set_title("SubX {var} vs. {comparison_var} {var}")
#             ax.set_xlabel("{comparison_var} {comparison_type}")
#             ax.set_ylabel("SubX {comparison_type}")
#             ax.set(xlim=(0,max_axis_val))
#             ax.set(ylim=(0,max_axis_val))
#             plt.plot(x_line,y_line,color='r')
#             ax.text(v1, v2,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=9) #add text
#             plt.savefig(f'{image_dir}/{var}_MSE.tif', dpi=300)
#             plt.close()
#         try:
#             plt.imread(f'{image_dir}/ETo_heatmap_1.tif')
#         except FileNotFoundError:
#             #Plot heatmap 1
#             heatmap, xedges, yedges = np.histogram2d(x, y, bins=100)
#             extent = [xedges[0], yedges[-1], xedges[0], yedges[-1]]
#             plt.imshow(heatmap.T, extent=extent, origin='lower')
#             plt.text(v3, v4,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=14, color = 'red')
#             plt.xlabel('{comparison_var} {var} {comparison_type}')
#             plt.ylabel('SubX {var} {comparison_type}')
#             plt.savefig(f'{image_dir}/{var}_heatmap_1.tif',dpi=300)
#             plt.close()

             
#         try:
#             plt.imread(f'{image_dir}/{var}_MSE_mean.tif')
#         except FileNotFoundError:
#             ax = sns.scatterplot(x = mean_x, y = mean_y,s=0.25)
#             ax.set_title("SubX {var} vs. {comparison_var} {var}")
#             ax.set_xlabel("{comparison_var} {comparison_type}")
#             ax.set_ylabel("SubX {comparison_type}")
#             ax.set(xlim=(0,max_axis_val))
#             ax.set(ylim=(0,max_axis_val))
#             plt.plot(x_line,y_line,color='r')
#             ax.text(v1, v2,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=9) #add text
#             plt.savefig(f'{image_dir}/{var}_MSE_mean.tif', dpi=300)
#             plt.close()
            
             
#         try:
#             plt.imread(f'{image_dir}/ETo_heatmap_mean.tif')
#         except FileNotFoundError:
#         #Plot heatmap
#             heatmap, xedges, yedges = np.histogram2d(mean_x, mean_y, bins=100)
#             extent = [xedges[0], yedges[-1], xedges[0], yedges[-1]]
#             plt.imshow(heatmap.T, extent=extent, origin='lower')
#             plt.text(v3, v5,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=14, color = 'red')
#             plt.xlabel('{comparison_var} {var} {comparison_type}')
#             plt.ylabel('SubX {var} {comparison_type}')
#             plt.savefig(f'{image_dir}/{var}_heatmap_1.tif',dpi=300)
#             plt.close()
   
            
# #Plot ETo
# plot_heatmaps_MSE(var='ETo',comparison_var='gridMET',comparison_type = 'mm/d',max_axis_val =30)


# def flatten_RZSM(date_list):
#     #Open SMERGE and SubX file
    
#     ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().values.ravel())    
#     file_length = len(date_list)
    
#     final_length = ravel_length * file_length
    
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_ETo = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open SMERGE file and ravel, and append to empty array
#         smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].values.ravel()
#         #append smerge to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_ETo[count_value:ravel_end,0] = smerge_file[:]
#         correlation_ETo[count_value:ravel_end,1] = sub_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
#     return(correlation_ETo)

# '''Find only the realization mean and plot scatterplot and heatmap'''

# def average_realization_RZSM(date_list):
#     #Open gridMET and SubX file
    
#     ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().mean(dim='model').values.ravel())    
#     file_length = len(date_list)
    
#     final_length = ravel_length * file_length
    
#     #Create an array with the length final_length and with 2 columns for gridMET and SubX
#     correlation_RZSM_realization_mean = np.empty((final_length, 2))

#     count_value = 0 #Keep a counter for index
#     ravel_end = ravel_length
#     for date_ in date_list:
        
#         #Open gridMET file and ravel, and append to empty array
#         #Open SMERGE file and ravel, and append to empty array
#         smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().mean(dim='model').values.ravel()
#         #Open SubX ETo file. Remove 1st dimension becuase it's useless
#         sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].mean(dim='model').values.ravel()
#         #append gridMET to first column (0 index)
#         #append SubX to 2nd columnn (1st index)
#         correlation_RZSM_realization_mean[count_value:ravel_end,0] = smerge_file[:]
#         correlation_RZSM_realization_mean[count_value:ravel_end,1] = sub_file[:]
        
#         #Update counter to append properly
#         count_value += ravel_length
#         ravel_end += ravel_length
        
#     return(correlation_RZSM_realization_mean)

        

# def SM_plot_setup(date_list):
    
#     correlation_SM = flatten_RZSM(date_list)    
#     correlation_RZSM_realization_mean = (average_realization_RZSM(date_list))

#     #Find the areas where there are only values in the gridMET dataset because
#     #gridMET is already restricted to CONUS.
#     #make inf values np.nan
#     correlation_SM[correlation_SM > 1e308] = np.nan
#     #Make all negative values equal to 0 because Ref. ET cannot be negative
#     correlation_SM[correlation_SM < 0] = 0
#     #remove gridMET grid cells with np.nan because of outside CONUS bbox
#     corr_1 = correlation_SM[~np.isnan(correlation_SM[:,0])]
#     #remove SubX values with np.nan because of differences in grid coordinates/missing values in a few CONUS grid cells
#     corr_2 = corr_1[~np.isnan(corr_1[:,1])]

   
#     #Same with mean
#     correlation_RZSM_realization_mean[correlation_RZSM_realization_mean > 1e308] = np.nan
#     #Make all negative values equal to 0 because Ref. ET cannot be negative
#     correlation_RZSM_realization_mean[correlation_RZSM_realization_mean < 0] = 0
    
#     no_nan_mean = correlation_RZSM_realization_mean[(~np.isnan(correlation_RZSM_realization_mean[:,0]) & (~np.isnan(correlation_RZSM_realization_mean[:,1]))) ]
        
    
    
#     # Scatterplot and Correlations
#     # Data
    
#     #GridMET is x
#     x = corr_2[:,0]
#     #SubX is y
#     y = corr_2[:,1]
        
#     mean_x = no_nan_mean[:,0]
#     mean_y = no_nan_mean[:,1]
    
#     return(x,y,mean_x,mean_y)



# #Plot RZSM
# plot_heatmaps_MSE(var='RZSM',comparison_var='SMERGE',comparison_type = 'm3/m3',max_axis_val =1.0)


# try:
#     plt.imread(f'{image_dir}/RZSM_MSE.tif')
#     plt.imread(f'{image_dir}/RZSM_heatmap.tif')
#     plt.imread(f'{image_dir}/RZSM_MSE_mean.tif')
#     plt.imread(f'{image_dir}/RZSM_heatmap_mean.tif')
    
# except FileNotFoundError:
#     print('Finding MSE and heatmap for RZSM.')    
#     a,b,mean_a,mean_b = SM_plot_setup(date_list)
    
    
#     #Plot RZSM MSE
#     max_RZSM_val = 1.0
    
#     x_line = np.linspace(0,max_RZSM_val,100)
#     y_line = np.linspace(0,max_RZSM_val,100)
    
#     try:
#         plt.imread(f'{image_dir}/RZSM_MSE.tif')
#     except FileNotFoundError:
#         ax = sns.scatterplot(x = a, y = b,s=0.25)
#         ax.set_title("SubX RZSM vs. SMERGE RZSM")
#         ax.set_xlabel("SMERGE (m3/m3)")
#         ax.set_ylabel("SubX (m3/m3)")
#         ax.set(xlim=(0,max_RZSM_val))
#         ax.set(ylim=(0,max_RZSM_val))
#         plt.plot(x_line,y_line,color='r')
#         ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(a,b),5)}", fontsize=9) #add text
#         plt.savefig(f'{image_dir}/RZSM_MSE.tif', dpi=300)
#         plt.close()
    
#     try:
#         plt.imread(f'{image_dir}/RZSM_heatmap.tif')
#     except FileNotFoundError:
#         #Plot heatmap
#         heatmap, xedges, yedges = np.histogram2d(a, b, bins=100)
#         extent = [0, yedges[-1], 0, yedges[-1]]
#         plt.imshow(heatmap.T, extent=extent, origin='lower')
#         plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(a,b),5)}", fontsize=14, color = 'red')
#         plt.xlabel('SMERGE RZSM (m3/m3)')
#         plt.ylabel('SubX RZSM (m3/m3)')
#         plt.savefig(f'{image_dir}/RZSM_heatmap.tif',dpi=300)
#         plt.close()
        
        
#     try:
#         plt.imread(f'{image_dir}/RZSM_MSE_mean.tif')
#     except FileNotFoundError:
#         ax = sns.scatterplot(x = mean_a, y = mean_b,s=0.25)
#         ax.set_title("SubX RZSM vs. SMERGE RZSM")
#         ax.set_xlabel("SMERGE (m3/m3)")
#         ax.set_ylabel("SubX (m3/m3)")
#         ax.set(xlim=(0,max_RZSM_val))
#         ax.set(ylim=(0,max_RZSM_val))
#         plt.plot(x_line,y_line,color='r')
#         ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(mean_a,mean_b),5)}", fontsize=9) #add text
#         plt.savefig(f'{image_dir}/RZSM_MSE_mean.tif', dpi=300)
#         plt.close()
        
         
#     try:
#         plt.imread(f'{image_dir}/RZSM_heatmap_mean.tif')
#     except FileNotFoundError:
#     #Plot heatmap
#         heatmap, xedges, yedges = np.histogram2d(mean_a, mean_b, bins=100)
#         extent = [0, yedges[-1], 0, yedges[-1]]
#         plt.imshow(heatmap.T, extent=extent, origin='lower')
#         plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(mean_a,mean_b),5)}", fontsize=14, color = 'red')
#         plt.xlabel('SMERGE RZSM (m3/m3)')
#         plt.ylabel('SubX RZSM (m3/m3)')
#         plt.savefig(f'{image_dir}/RZSM_heatmap_mean.tif',dpi=300)
#         plt.close()        
        
    
#%%POTENTIAL PLOTS THAT DIDN'T WORK OUT  --heatmap will work, plt.scatter will not

# #This will eat up all the memory memory
# # Plot
# plt.rcParams.update({'figure.figsize':(10,8), 'figure.dpi':100})
# plt.scatter(x, y, label=f'y1 Correlation = {np.round(np.corrcoef(x,y)[0,1], 4)}')
# plt.title('Scatterplot and Correlations')
# plt.legend()
# plt.show()



    
    
    
    
    
    
    
    
    

