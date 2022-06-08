#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Goal is to look at each SubX file for soil moisture and reference ET and look at the
scatterplot between observed (soil moisture - SMERGE) and observed (reference ET - gridMET METDATA)
and between SubX files of the same variable.

File list dates = actual dates values start at 1-11-1999 and end at 2-09-2016



@author: kdl
"""



#TODO : Add a try and except clause to for when the file has already been completed


import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from multiprocessing import Pool
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error

#open a single file for SubX reference ET and a single file from gridMET 
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = int('10')

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

subX_dir = f'{dir1}/Data/SubX/GMAO'
subX_SM_out_dir = f'{dir1}/Data/SubX/GMAO/SM_converted_m3_m3' #conversion

gridMET_dir = f'{dir1}/Data/gridMET' #reference evapotranspiration
output_ETo_dir = f'{gridMET_dir}/ETo_SubX_values' #Refernce ET output directory
smerge_in_dir = f'{dir1}/Data/SMERGE_SM/Raw_data' #raw files that have been boxed in CONUS
SM_SubX_out_dir = f'{dir1}/Data/SMERGE_SM/SM_SubX_values' #Smerge values overlayed on SubX grid
eddi_dir = f'{dir1}/Data/EDDI/convert_2_nc'
eddi_subX_dir = f'{dir1}/Data/EDDI/EDDI_SubX_values'

image_dir = f'{dir1}/Outputs/MSE_plots'

#Make new directories
os.system(f'mkdir -p {image_dir}')
os.system(f'mkdir {SM_SubX_out_dir}')
os.system(f'mkdir {eddi_subX_dir}')
os.system(f'mkdir {output_ETo_dir}')

os.chdir(subX_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob('mrso*_*-*.nc'))
date_list = [i[-13:-3] for i in date_list]

# _date = date_list[0]

#%%
def EDDI_SubX_creation(_date):    
    try:
        xr.open_dataset(f'{eddi_subX_dir}/EDDI_SubX_{_date}.nc')
        print(f'Already completed date {_date}. Saved in {eddi_subX_dir}')
    except FileNotFoundError:
        
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        eddi_file = xr.open_dataset(f'{eddi_dir}/eddi_merged.nc4').astype('float64')
        eddi_file = eddi_file.sel(time=slice("1999-01-10","2016-02-09"))

        print(f'Working on initialized day {_date} to retrieve EDDI values from SubX models, leads, & coordinates and saving data into {eddi_subX_dir}.')
        var_file = xr.open_dataset(f'{subX_dir}/EDDI_{_date}.nc4')
        
        #Remove dimension S[1] because it has not data (it's a GMAO GEOS specific issue)
        var_file_arr = var_file.EDDI[:,:,:,:,:].values
        
        
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
        eddi_out = (np.zeros_like(var_file_arr).squeeze()).astype('float64')
     
        for model in range(var_file.EDDI.shape[1]):
            for i_lead in range(var_file.EDDI.shape[2]):
                date_val = pd.to_datetime(var_file.S.values[0]) + dt.timedelta(days=i_lead)
            
                for i_Y in range(var_file.EDDI.shape[3]):
                    for i_X in range(var_file.EDDI.shape[4]):
                        
                        eddi_out[0,model, i_lead, i_Y, i_X] = \
                            eddi_file.EDDI.sel(time = date_val).isel(X = i_X, Y = i_Y).values
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_SubX_value = (['S','model','lead','Y','X'], eddi_out[:,:,:,:,:]),
            ),
            coords = dict(
                S = var_file.S.values,
                X = var_file.X.values,
                Y = var_file.Y.values,
                lead = var_file.lead.values,
                model = var_file.model.values,
        
            ),
            attrs = dict(
                Description = 'Evaporative Demand Drought Index values on the exact same date and grid \
                cell as SubX data'),
        )                    
            
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{eddi_subX_dir}/ETo_SubX_{_date}.nc4', mode ='w')
        #To find the correct date
        
        print(f'Completed ETo_SubX_{_date}')
print(f'Working from directory {subX_dir}')


#%% ETo gridMET
def ETo_gridMET_SubX_creation(_date) -> float:    
    try:
        xr.open_dataset(f'{output_ETo_dir}/ETo_SubX_{_date}.nc')
        print(f'Already completed date {_date}. Saved in {output_ETo_dir}')
    except FileNotFoundError:
        
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        gridMET_file = xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc').astype('float64')
        gridMET_file = gridMET_file.sel(day=slice("1999-01-10","2016-02-09"))

        print(f'Working on initialized day {_date} to find gridMET values from SubX models, leads, & coordinates and saving data into {output_ETo_dir}.')
        sub_file = xr.open_dataset(f'ETo_{_date}.nc4')
        
        #Remove dimension S[1] because it has not data (it's a GMAO GEOS specific issue)
        sub_file_arr = sub_file.ETo[:,:,:,:,:].values
        
        
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
        gridMET_out = (np.zeros_like(sub_file_arr).squeeze()).astype('float64')
     
        for model in range(sub_file.ETo.shape[1]):
            for i_lead in range(sub_file.ETo.shape[2]):
                date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead)
            
                for i_Y in range(sub_file.ETo.shape[3]):
                    for i_X in range(sub_file.ETo.shape[4]):
                        
                        gridMET_out[0,model, i_lead, i_Y, i_X] = \
                            gridMET_file.ETo_gridmet.sel(day = date_val).isel(lon = i_X, lat = i_Y).values
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_SubX_value = (['S','model','lead','Y','X'], gridMET_out[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = sub_file.lead.values,
                model = sub_file.model.values,
        
            ),
            attrs = dict(
                Description = 'gridMET Reference ETo values on the exact same date and grid \
                cell as SubX data'),
        )                    
            
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_ETo_dir}/ETo_SubX_{_date}.nc', mode ='w')
        #To find the correct date
        
        print(f'Completed ETo_SubX_{_date}')

#%%#SMERGE Soil Moisture
#SubX.dat depth file
depth_file = np.loadtxt(f'{dir1}/Data/SubX/SubX_sm_grid_conversion/prof_depth_1x1.dat')

'''As opposed to just shifting the coordinates, we need to take the average of the 
the longitude that is closest to the SubX longitude value'''

#Change Longitude to 360 degree notation
for row in range(len(depth_file)):
    if depth_file[row,0] < 0:
        depth_file[row,0] = 360 + depth_file[row,0]
            
#Latitude is already on a 180 degree going from -89.5 to 89.5
    
def convert_SubX_to_m3_m3(_date):
    try:
       xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{_date}.nc4')
       print(f'Already completed date {_date}. Saved in {subX_SM_out_dir}')
       
    except FileNotFoundError:
        
        print(f'Working on initialized day {_date} to convert SubX soil moisture into m3/m3 and saving into {subX_SM_out_dir}.')
            
        '''To convert from SUBX kg/m2 (area), we need to divide each value based 
    on latitude and longitude by the SubX.dat file which contains the depth
    of soil column in meters (which we convert to mm).'''
                
        global depth_file
        
        SubX_file = xr.open_dataset(f'{subX_dir}/mrso_GMAO_{_date}.nc')
        
        out_sm_SubX = xr.zeros_like(SubX_file)
        

        for model in range(SubX_file.M.shape[0]):
            for lead in range(SubX_file.L.shape[0]):
                for Y in range(SubX_file.Y.shape[0]):
                    for X in range(SubX_file.X.shape[0]):
                        
                        #Find what the depth Y and X coordinates would be on either end by 1/2 degree                    
                        Y_plus_value = SubX_file.Y[Y].values + 0.5 
                        Y_minus_value = SubX_file.Y[Y].values - 0.5
                        
                        X_plus_value = SubX_file.X[X].values + 0.5 
                        X_minus_value = SubX_file.X[X].values - 0.5

                        val_plus_index = np.where((depth_file[:,0] == X_plus_value) & (depth_file[:,1] == Y_plus_value))
                        val_plus = depth_file[val_plus_index][0,-1]

                        val_minus_index = np.where((depth_file[:,0] == X_minus_value) & (depth_file[:,1] == Y_minus_value))
                        val_minus = depth_file[val_minus_index][0,-1]
                        
                       
                        #Dont work on empty grid cells
                        if np.count_nonzero(np.isnan(SubX_file.mrso[0,model, lead, Y, X].values)) == 1:
                            val_3 = np.nan
                        
                        #Avearge of lon or lat due to grid differences
                        #Only look at specific ranges of CONUS
 
                            '''If there are no elevation values for a specific grid cell and the grid cell
                        is supposed to be filled based on only soil moisture value, take the average
                        of the nearest 8 block of grid cells'''
                        
                        elif np.nansum(val_plus) == 0 and np.nansum(val_minus) == 0 and \
                            np.count_nonzero(np.isnan(SubX_file.mrso[0,model, lead, Y, X].values)) != 1:
                            #Find nearby elevation values
                            #Subtract 1 from latisubX_SM_out_dirtude
                            y_m1 = np.where((depth_file[:,0] == X_plus_value) & (depth_file[:,1] == Y_plus_value - 1))
                            y_m1 = depth_file[y_m1][0,-1]
                            #Add one to latitude
                            y_a1 = np.where((depth_file[:,0] == X_plus_value) & (depth_file[:,1] == Y_plus_value + 1))
                            y_a1 = depth_file[y_a1][0,-1]
                            
                            #Subtract 1 from longitude
                            x_m1 = np.where((depth_file[:,0] == X_plus_value -1) & (depth_file[:,1] == Y_plus_value))
                            x_m1 = depth_file[x_m1][0,-1]
                            #Add one to latitude
                            x_a1 = np.where((depth_file[:,0] == X_plus_value + 1) & (depth_file[:,1] == Y_plus_value))
                            x_a1 = depth_file[x_a1][0,-1]            
                            
                            #Subtract 1 from latitude; add 1 to latitude
                            yx_m1 = np.where((depth_file[:,0] == X_plus_value + 1) & (depth_file[:,1] == Y_plus_value - 1))
                            yx_m1 = depth_file[yx_m1][0,-1]
                            #Add one to latitude
                            yx_a1 = np.where((depth_file[:,0] == X_plus_value -1 ) & (depth_file[:,1] == Y_plus_value + 1))
                            yx_a1 = depth_file[yx_a1][0,-1]
                            
                            #Subtract 1 from longitude
                            xy_m1 = np.where((depth_file[:,0] == X_plus_value -1) & (depth_file[:,1] == Y_plus_value + 1))
                            xy_m1 = depth_file[xy_m1][0,-1]
                            #Add one to latitude
                            xy_a1 = np.where((depth_file[:,0] == X_plus_value + 1) & (depth_file[:,1] == Y_plus_value -1))
                            xy_a1 = depth_file[xy_a1][0,-1]            
                            
                            
                            total_boxes = 8
                            miss_data = np.count_nonzero(np.isnan((y_m1, y_a1, x_m1, x_a1, yx_m1, yx_a1, xy_m1, xy_a1)))
                            
                            dvisor = total_boxes - miss_data
                            
                            val_3 = np.nansum((y_m1, y_a1, x_m1, x_a1, yx_m1, yx_a1, xy_m1, xy_a1))/ dvisor
                            
                        elif np.nansum(val_plus) == 0:
                            val_3 = val_minus
                        elif np.nansum(val_minus) == 0:
                            val_3 = val_plus
                        else:
                            val_3 = np.add(val_plus,val_minus)/2

                        #Divide by lat/lon depth in mm (current values in data are in meters)
                        #multiply meters by 1000 to get mm
                        out_sm_SubX.mrso[0,model, lead, Y, X] =  \
                            SubX_file.mrso[0,model, lead, Y, X].values / ((val_3)*1000)
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
        data_vars = dict(
            SM_SubX_m3_m3_value = (['S','model','lead','Y','X'], out_sm_SubX.mrso[:,:,:,:,:]),
        ),
        coords = dict(
        S = out_sm_SubX.S.values,
        X = out_sm_SubX.X.values,
        Y = out_sm_SubX.Y.values,
        lead = out_sm_SubX.L.values,
        model = out_sm_SubX.M.values,
        ),
        attrs = dict(
        Description = 'SubX SM (kg/m2) converted to (m3/m3) based on profile depth.'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{subX_SM_out_dir}/SM_SubX_m3_m3_{_date}.nc4', mode ='w')
        print(f'Completed SM_m3_m3_{_date}.')
    

#%%
    
def SM_SMERGE_SubX_creation(_date):    
    try:
        xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{_date}.nc')
        print(f'Already completed date {_date}. Saved in {SM_SubX_out_dir}.')
        
    except FileNotFoundError:
            
        print(f'Working on initialized day {_date} to find SMERGE values from SubX models, leads, & coordinates and saving data into {SM_SubX_out_dir}.')
        sub_file = xr.open_dataset(f'mrso_GMAO_{_date}.nc')
        
        #Remove dimension S[1] because it has not data (it's a GMAO GEOS specific issue)
        sub_file_arr = sub_file.mrso[:,:,:,:,:].values
        
        SMERGE_file = xr.open_dataset(f'{smerge_in_dir}/smerge_sm_merged_remap.nc4')
        
        '''Same format as ETo. Find the correct dates, leads, models and just fill in'''
        smerge_out = (np.zeros_like(sub_file_arr).squeeze()).astype('float64')
     
        for model in range(sub_file.mrso.shape[1]):
            for i_lead in range(sub_file.mrso.shape[2]):
                date_val = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead)
            
                for i_Y in range(sub_file.mrso.shape[3]):
                    for i_X in range(sub_file.mrso.shape[4]):
                        
                        smerge_out[0,model, i_lead, i_Y, i_X] = \
                            SMERGE_file.RZSM.sel(time = date_val).isel(X = i_X, Y = i_Y).values
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                SMERGE_SubX_value = (['S','model','lead','Y','X'], smerge_out[:,:,:,:,:]),
            ),
            coords = dict(
                S = sub_file.S.values,
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                lead = sub_file.L.values,
                model = sub_file.M.values,
        
            ),
            attrs = dict(
                Description = 'SMERGE SM (m3/m3) values on the exact same date and grid \
                cell as SubX data'),
        )                    
            
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{SM_SubX_out_dir}/SM_SubX_{_date}.nc', mode ='w')
        #To find the correct date
        print(f'Completed SM_SubX_{_date}')        


#%%Run all four functions
if __name__ == '__main__':
    print('Converting SubX SM to m3/m3.')
    p = Pool(num_processors)
    p.map(convert_SubX_to_m3_m3, date_list)
    p.map(SM_SMERGE_SubX_creation, date_list)
    print('Working on EDDI.')
    p.map(ETo_gridMET_SubX_creation, date_list)
    print('Working on gridMET ETo')
    p.map(EDDI_SubX_creation, date_list)

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



    
    
    
    
    
    
    
    
    

