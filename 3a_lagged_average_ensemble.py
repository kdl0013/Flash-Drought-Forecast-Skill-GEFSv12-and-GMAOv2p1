#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:06:40 2022

@author: Kyle Lesinger

Testing how to make a lagged ensemble average away from Cheyenne cluster

Initial test:
    Create a new file that contains stacked dimensions that would contain a larger
    distribution for soil moisture percentile. Issue was >500GB memory is needed
    for the 6 dimensional dataset I was trying to produce.
    
New test:
    Just make a lagged average ensemble. This will create a smaller distribution
    for soil moisture percentiles, but should still have a plenty large enough 
    distribution to pick from.
    
Inputs:
    Individual SubX files for each variable. 
    
Output:
    Lagged average ensemble for each model and variable.

"""

import xarray as xr
import numpy as np
import os
from glob import glob
import pandas as pd
import datetime
from collections import OrderedDict
import time 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error

model = 'GMAO'

home = 'main_dir'
# home='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

variables = ['SM_SubX_m3_m3', 'ETo']

home_dir = f'{home}/Data/SubX/{model}'
sm_dir = f'{home_dir}/SM_converted_m3_m3'
smerge_dir = f'{home}/Data/SMERGE_SM/Raw_data'
image_dir = f'{home}/Outputs/MSE_plots'
gridMET_dir = f'{home}/Data/gridMET'


os.chdir(home_dir)
os.system('mkdir lagged_ensemble_by_variable')

#Must manually change certain values to the actual variable shorthand because
#xarray needs 
#%%
for var in variables:

    if var == 'SM_SubX_m3_m3':
        dir_ = sm_dir
        var_name = 'SM_SubX_m3_m3_value'
    elif var == 'ETo':
        dir_ = home_dir
        var_name = var

    
    print(f'Starting variable {var}. Sleeping 5 seconds. This will take about 15 minutes')
    time.sleep(5)
    try:
        #Do not complete script if file calculation has already been conducted
        var_completed_file = xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/{var}_lagged_average.nc4')
    except FileNotFoundError:
            
        #get info for future dimension creation later
        def inspect_files(var, date_start_val=-14, date_end_val=-4):
            '''Pre-processing, open the first file and inspect contents, save date list
            Inputs: 
                var = variable, date_start_val is the index location of date in file name
            Outputs: 
                  list of dates in directory (date_list) and a single file to inspect (varFile) 
                  and a count of files for that variable in directory. 
            '''
            
    
            
            count=0
            date_list = []
            for file in sorted(glob(f'{dir_}/{var}_*.nc4')):
                if count == 0:
                    #print(file)
                    varFile = xr.open_dataset(file)
                    date_list.append(file[date_start_val:date_end_val])
                else:
                    date_list.append(file[date_start_val:date_end_val])
                    
                count+=1
                
            # get the dates for all SubX models (this is all the possible dates)
            start_date, end_date = pd.to_datetime((date_list[0])), pd.to_datetime((date_list[-1]))
            
            #dates
            daily_dates = [start_date + datetime.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 46)] #Add final initialization days dates
            
            return(date_list, varFile, count, daily_dates, varFile)  
        
        
                
        date_list, varFile, count_files, daily_dates, varFileTEST = inspect_files(var=var)
        #Convert date_list to a daily time series
        
        varFile.info()
        
        #Get the difference between ordinal dates
        date_ordinal = []
        for i in date_list:
            date_ordinal.append(pd.to_datetime(i).toordinal())
        
        diff_days = np.diff(np.array(date_ordinal))
        np.unique(diff_days) #Difference between days, 1 file has a 35 day difference at index
        np.where(diff_days >6) #at index 916
        
        #expand dimensions function for xarray
        def expand_dimensions(data, fill_value=np.nan, **new_coords):
            """Expand (or add if it doesn't yet exist) the data array to fill in new
            coordinates across multiple dimensions.
        
            If a dimension doesn't exist in the dataarray yet, then the result will be
            `data`, broadcasted across this dimension.
        
            >>> da = xr.DataArray([1, 2, 3], dims="a", coords=[[0, 1, 2]])
            >>> expand_dimensions(da, b=[1, 2, 3, 4, 5])
            <xarray.DataArray (a: 3, b: 5)>
            array([[ 1.,  1.,  1.,  1.,  1.],
                   [ 2.,  2.,  2.,  2.,  2.],
                   [ 3.,  3.,  3.,  3.,  3.]])
            Coordinates:
              * a        (a) int64 0 1 2
              * b        (b) int64 1 2 3 4 5
        
            Or, if `dim` is already a dimension in `data`, then any new coordinate
            values in `new_coords` that are not yet in `data[dim]` will be added,
            and the values corresponding to those new coordinates will be `fill_value`.
        
            >>> da = xr.DataArray([1, 2, 3], dims="a", coords=[[0, 1, 2]])
            >>> expand_dimensions(da, a=[1, 2, 3, 4, 5])
            <xarray.DataArray (a: 6)>
            array([ 1.,  2.,  3.,  0.,  0.,  0.])
            Coordinates:
              * a        (a) int64 0 1 2 3 4 5
        
            Args:
                data (xarray.DataArray):
                    Data that needs dimensions expanded.
                fill_value (scalar, xarray.DataArray, optional):
                    If expanding new coords this is the value of the new datum.
                    Defaults to `np.nan`.
                **new_coords (list[int | str]):
                    The keywords are arbitrary dimensions and the values are
                    coordinates of those dimensions that the data will include after it
                    has been expanded.
            Returns:
                xarray.DataArray:
                    Data that had its dimensions expanded to include the new
                    coordinates.
            """
            ordered_coord_dict = OrderedDict(new_coords)
            shape_da = xr.DataArray(
                np.zeros(list(map(len, ordered_coord_dict.values()))),
                coords=ordered_coord_dict,
                dims=ordered_coord_dict.keys())
            expanded_data = xr.broadcast(data, shape_da)[0].fillna(fill_value)
            return expanded_data
        

        #Get info about date list and the actual number of days
        number_days = pd.to_datetime(date_list[-1]).toordinal() - pd.to_datetime(date_list[0]).toordinal()
        print()
        print()
        print()
        print()
        
        print('######### Number of Days in Final Time Series (Adding 45 days to end file)')
        print(f'Total number of days in entire time series (plus lead days) is {number_days+varFile.lead.shape[0] + 1}. Day 0 has no values though.')
        #%%Find the date with values - completed, S[0] dimension has the data
        '''Some files have 2 S dimensions which means it has 2 initialization days.
        Upon further inspection, it appears that some files have no values in one of the 
        dimensions. GMAO GEOS model specifically.
        
        The first observed date (1999-01-10.nc) has 2 dimensions. The original download
        says that it was initialized on 1999-01-10, but the actual S value date is 
        1999-01-11.  For this first file, the S dimension with values is S[0] {1999-01-11} 
        and not S[1] {1999-01-10}. 
        
        After inspection of all files, dwsrf file name is actually 1 day behind the actual 
        initialization date. In other words, the S[0] dimension value (the date) actual is 1 day 
        ahead of the file name.
        
        After further thinking, it appears that the lead of 0 is actually a day's lead (24 hours)
        
        '''
        print()
        print()
        print()
        print()
        if var == 'SM_SubX_m3_m3':
            dim_0 = 0
            dim_1 = 1
            for file in sorted(glob(f'{dir_}/{var}_*.nc4'))[0:9]:
                print(f'File date says:  {file[-14:-4]}')
                file_o = xr.open_dataset(file)
                print(f'Number of missing values : {str(np.count_nonzero(np.isnan(file_o[var_name][dim_0,:,:,:,:].values)))}' + f'   Values are in S[{dim_0}] dim   Day {str(file_o.S.values[dim_0])[0:10]}') 
                print(f'Number of missing values : {str(np.count_nonzero(np.isnan(file_o[var_name][dim_1,:,:,:,:].values)))}' + f'    Values are in S[{dim_1}] dim   Day {str(file_o.S.values[dim_1])[0:10]}') 
            
        #You can visually see that the first dimension has no np.nan, but S[1] dimension does
        #%%
        '''Create a new dataset to store final output'''
        #Step 2, open a file and get the dates for each S dimension as well as create
        #more dates based on lead time. 
        
        #Stack a single file grid and make a copy to create a new empty dataset
        varFile = varFile.stack(grid = ['Y','X'])
        var_F_unstack = varFile.unstack('grid')
        
        #Continue settin up empty dataset to append lagged average to.
        varFile = (varFile[var_name][0,:,0,:]).to_dataset()
        varFile = varFile.drop(labels = ['lead','model','S'])
        
        #Add 3rd dimension for time series dates
        varFile = expand_dimensions(varFile, time_series_dates = pd.to_datetime(daily_dates))

        varFile[var_name][:,:,:] = (0)
        varFile= varFile.astype('float64')
        
        # varFile = varFile.unstack('grid')
        var_TEST = varFile.copy()
        #Output shape, use var_TEST as the output for the function
        var_TEST[var_name].shape
        #convert to numpy array
        var_TEST = var_TEST.unstack('grid')
        var_TEST = np.array(var_TEST.to_array())
        var_TEST = var_TEST.squeeze()
        #For finding number of values to take correct lagged average
        num_values_TEST = np.zeros_like(var_TEST)

#%%        
        '''LAGGED AVERAGE ENSEMBLE FILE CREATION
        
        Input: 
            File names
        Output:
            [4 (realization) x 1593 (grid) x time_series length]
            var_TEST array which contains the average value for each day's file
            when by filling in the ordinal dates (file_date - beginning_file_date)
            
            This accounts for the different models and files which each have
            different days associated with each file.
        '''
        
        for file in sorted(glob(f'{dir_}/{var}_*.nc4')):
            #Convert to numpy array, but keep date information
            def file_setup(file):
                #Read file
                open_f = xr.open_dataset(file)
                #Get date information
                beginning_start_date_ordinal = pd.to_datetime(list(sorted(glob('*.nc4')))[0][-14:-4], format= '%Y/%m/%d').toordinal()
                beginning_start_date = pd.to_datetime(list(sorted(glob('*.nc4')))[0][-14:-4])
                beginining_index = date_list.index((list(sorted(glob('*.nc4')))[0][-14:-4]))
                
                open_f_ordinal = pd.to_datetime(file[-14:-4], format= '%Y/%m/%d').toordinal()
                #To keep up with file number within dimensions
                file_num = list(sorted(glob(f'{dir_}/*.nc4'))).index(file) #Number of files. Each file in its own dimension (file_size_list)
                
                #Get date info about opened, intialized file
                open_f_lead_0_date = str(pd.to_datetime(open_f.S.values[0], format= '%Y/%m/%d'))[0:10] #Character string of lead 0 day (day after initialization)
                open_f_lead_0_minus_1_date = str(pd.to_datetime(open_f.S.values[1], format= '%Y/%m/%d'))[0:10]
                
                open_f_lead_0_ordinal = pd.to_datetime(open_f_lead_0_date, format= '%Y/%m/%d').toordinal()
                open_f_lead_0_minus_1_ordinal = pd.to_datetime(open_f_lead_0_minus_1_date, format= '%Y/%m/%d').toordinal()
                
                open_f_lead_0_minus_1_index = file_num  #Same as file n
                
                #Convert to numpy array
                open_f_arr = open_f.to_array()
                open_f_arr = np.array(open_f_arr)
                #Drop the first dimension (contains nothing)
                open_f_arr = open_f_arr.squeeze()
                #Keep only the first axis values (actually contains values)
                open_f_arr = open_f_arr[0,:,:,:]
                #Reshape
                # open_f_arr = np.reshape(open_f_arr,(open_f_arr.shape[0],open_f_arr.shape[2],open_f_arr.shape[1]))
                #Now we can maybe add Numba to fill in the missing values
                vector = np.vectorize(float)
                open_f_arr = vector(open_f_arr)
                
                return(open_f_arr, int(open_f_lead_0_ordinal), \
                        int(beginning_start_date_ordinal), int(file_num))
            
            open_f_arr, open_f_lead_0_ordinal, beginning_start_date_ordinal, file_num = file_setup(file)
        
            
            #Take the average of each grid cell with it's recursive date and input
            #into var_TEST as the output
            # @vectorize([int64(int64,int64)])
            def lagged_ensemble_preprocess(open_f_lead_0_ordinal, beginning_start_date_ordinal):
                global open_f_arr
                global var_TEST
                global num_values_TEST
                
                for mod in range(open_f_arr.shape[0]):
                    for lead in range(open_f_arr.shape[1]):
                        lead_final_ordinal = open_f_lead_0_ordinal + lead
                    # date_count = 0 #For understanding the starting date of the file
                        for Yx in range(open_f_arr.shape[2]): 
                            for Xx in range(open_f_arr.shape[3]):
                                var_TEST[mod,(lead_final_ordinal - beginning_start_date_ordinal),Yx,Xx] = \
                                    var_TEST[mod,(lead_final_ordinal - beginning_start_date_ordinal),Yx,Xx] + \
                                        open_f_arr[mod,lead, Yx,Xx]
                                
                                #Checks the number of values for each grid cell 
                                #so we can take the average
                                if open_f_arr[mod,lead, Yx,Xx] > 0.0:
                                    num_values_TEST[mod,(lead_final_ordinal - beginning_start_date_ordinal),Yx,Xx] += 1

            #Run function
            lagged_ensemble_preprocess(open_f_lead_0_ordinal, beginning_start_date_ordinal)
        
        
        print(f'Completed lagged ensemble average for {var}')
        print()
        print()
        print()
        print()

# Save file to netcdf for later processing
        '''Need to merge soil moisture to get into right format for mask of CONUS
        Apply to only the mrso '''
        #Take the average to complete lagged ensemble average
        lagged_average_xr = var_TEST / num_values_TEST


        if var == 'mrso':
            desc = 'Vertically integrated Soil Moisture'
        elif var == 'tdps':
            desc = '2-meter Dew point temperature'
        elif var == 'tas':
            desc = '2-meter Air temperature'
        elif var == 'pr':
            desc = 'Total precipitation'
        elif var == 'uas':
            desc = '10-meter U component of wind'
        elif var == 'vas':
            desc = '10-meter V component of wind'
        elif var == 'cape':
            desc = 'Convective Available Potential Energy'
        elif var == 'dswrf':
            desc = 'Downwelling shortwave radiation'
        elif var == 'huss':
            desc = 'Specific humidity at 850hPa'
        elif var == 'SM_SubX_m3_m3':
            desc = 'SM_SubX_m3_m3'
        elif var ==  'ETo':
            desc = 'ETo_SubX_mm_d'


        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                READ_DESCRIPTION_ATTRS = (['realization','time_series_dates','Y','X'], lagged_average_xr),
            ),
            coords = dict(
                X = var_F_unstack.X.values,
                Y = var_F_unstack.Y.values,
                time_series_dates = daily_dates,
                realization = var_F_unstack.model.values,
            ),
            attrs = dict(
                Description = desc),
        )
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'lagged_ensemble_by_variable/{var}_lagged_average.nc4', mode ='w')
        
#%%Scatterplot and heatmap lagged average

def flatten_RZSM():
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    ravel_length = len(xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/SM_SubX_m3_m3_lagged_average.nc4').to_array().values.ravel())    
    ravel_mean = len(xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/SM_SubX_m3_m3_lagged_average.nc4').to_array().mean(dim='realization').values.ravel())    
    
    correlation_RZSM = np.empty((ravel_length, 2))
    correlation_RZSM_mean = np.empty((ravel_mean, 2))

    #Open subX lagged average file and smerge file
    sm_SubX_file = xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/SM_SubX_m3_m3_lagged_average.nc4')
    sm_SMERGE_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
    
    #First make the mean file
    SubX_sm_mean = sm_SubX_file.to_array().mean(dim='realization').values
    
    #Add raveled files to correlation_RZSM_mean
    #column 0 is SMERGE
    #column 1 is SubX
    correlation_RZSM_mean[:,0] = sm_SMERGE_file.RZSM.values.ravel()
    correlation_RZSM_mean[:,1] = SubX_sm_mean.ravel()
    
    no_nan_mean = correlation_RZSM_mean[(~np.isnan(correlation_RZSM_mean[:,0]) & (~np.isnan(correlation_RZSM_mean[:,1])))]
    mean_x,mean_y = no_nan_mean[:,0], no_nan_mean[:,1]
    
    #Now find the correlation between the lagged average file with all models and 
    #Smerge - add to correlation_RZSM
    ravel_init = 0
    ravel_len = int(ravel_length/sm_SubX_file.realization.shape[0])
    ravel_start = ravel_len
    for model in range(sm_SubX_file.realization.shape[0]):
        #ravel the first model, add to correlation_RZSM index column 1
        #index column 0 will be smerge
        
        
        correlation_RZSM[ravel_init:ravel_start,0] = sm_SMERGE_file.RZSM.values.ravel()
        correlation_RZSM[ravel_init:ravel_start,1] = sm_SubX_file.READ_DESCRIPTION_ATTRS[model,:,:,:].values.ravel()
        
        ravel_init += ravel_len
        ravel_start += ravel_len
    
    no_nan = correlation_RZSM[(~np.isnan(correlation_RZSM[:,0]) & (~np.isnan(correlation_RZSM[:,1])))]
    x,y = no_nan[:,0], no_nan[:,1]
    
    return(x,y,mean_x,mean_y)
    
#Plot heatmaps and scatterplots

try:
    plt.imread(f'{image_dir}/RZSM_MSE_lagged_average.tif')
    plt.imread(f'{image_dir}/RZSM_heatmap_lagged_average.tif')
    plt.imread(f'{image_dir}/RZSM_MSE_lagged_average_mean.tif')
    plt.imread(f'{image_dir}/RZSM_heatmap_lagged_average_mean.tif')
    
except FileNotFoundError:
    print('Finding MSE and heatmap for RZSM lagged average ensemble.')    
    x,y,mean_x,mean_y = flatten_RZSM()    
    
    #Plot RZSM MSE
    max_RZSM_val = 1.0
    
    x_line = np.linspace(0,max_RZSM_val,100)
    y_line = np.linspace(0,max_RZSM_val,100)
    
    try:
        plt.imread(f'{image_dir}/RZSM_MSE_lagged_average.tif')
    except FileNotFoundError:
        ax = sns.scatterplot(x = x, y = y,s=0.25)
        ax.set_title("SubX RZSM vs. SMERGE RZSM")
        ax.set_xlabel("SMERGE (m3/m3)")
        ax.set_ylabel("SubX (m3/m3)")
        ax.set(xlim=(0,max_RZSM_val))
        ax.set(ylim=(0,max_RZSM_val))
        plt.plot(x_line,y_line,color='r')
        ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=9) #add text
        plt.savefig(f'{image_dir}/RZSM_MSE_lagged_average.tif', dpi=300)
        plt.close()
    
    try:
        plt.imread(f'{image_dir}/RZSM_heatmap_lagged_average.tif')
    except FileNotFoundError:
        #Plot heatmap
        heatmap, xedges, yedges = np.histogram2d(x,y, bins=100)
        extent = [0, yedges[-1], 0, yedges[-1]]
        plt.imshow(heatmap.T, extent=extent, origin='lower')
        plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=14, color = 'red')
        plt.xlabel('SMERGE RZSM (m3/m3)')
        plt.ylabel('SubX RZSM (m3/m3)')
        plt.savefig(f'{image_dir}/RZSM_heatmap_lagged_average.tif',dpi=300)
        plt.close()
        
        
    try:
        plt.imread(f'{image_dir}/RZSM_MSE_lagged_average_mean.tif')
    except FileNotFoundError:
        ax = sns.scatterplot(x = mean_x, y = mean_y,s=0.25)
        ax.set_title("SubX RZSM vs. SMERGE RZSM")
        ax.set_xlabel("SMERGE (m3/m3)")
        ax.set_ylabel("SubX (m3/m3)")
        ax.set(xlim=(0,max_RZSM_val))
        ax.set(ylim=(0,max_RZSM_val))
        plt.plot(x_line,y_line,color='r')
        ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=9) #add text
        plt.savefig(f'{image_dir}/RZSM_MSE_lagged_average_mean.tif', dpi=300)
        plt.close()
        
         
    try:
        plt.imread(f'{image_dir}/RZSM_heatmap_lagged_average_mean.tif')
    except FileNotFoundError:
    #Plot heatmap
        heatmap, xedges, yedges = np.histogram2d(mean_x, mean_y, bins=100)
        extent = [0, yedges[-1], 0, yedges[-1]]
        plt.imshow(heatmap.T, extent=extent, origin='lower')
        plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=14, color = 'red')
        plt.xlabel('SMERGE RZSM (m3/m3)')
        plt.ylabel('SubX RZSM (m3/m3)')
        plt.savefig(f'{image_dir}/RZSM_heatmap_lagged_average_mean.tif',dpi=300)
        plt.close()        
        
        
#%%
def flatten_ETO():
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    ravel_length = len(xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/ETo_lagged_average.nc4').to_array().values.ravel())    
    ravel_mean = len(xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/ETo_lagged_average.nc4').to_array().mean(dim='realization').values.ravel())    
    
    correlation_RZSM = np.empty((ravel_length, 2))
    correlation_RZSM_mean = np.empty((ravel_mean, 2))

    #Open subX lagged average file and smerge file
    sm_SubX_file = xr.open_dataset(f'{home_dir}/lagged_ensemble_by_variable/ETo_lagged_average.nc4')
    sm_SMERGE_file = xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
    
    #First make the mean file
    SubX_sm_mean = sm_SubX_file.to_array().mean(dim='realization').values
    
    #Add raveled files to correlation_RZSM_mean
    #column 0 is SMERGE
    #column 1 is SubX
    correlation_RZSM_mean[:,0] = sm_SMERGE_file.ETo_gridmet.values.ravel()
    correlation_RZSM_mean[:,1] = SubX_sm_mean.ravel()
    
    no_nan_mean = correlation_RZSM_mean[(~np.isnan(correlation_RZSM_mean[:,0]) & (~np.isnan(correlation_RZSM_mean[:,1])))]
    mean_x,mean_y = no_nan_mean[:,0], no_nan_mean[:,1]
    
    #Now find the correlation between the lagged average file with all models and 
    #Smerge - add to correlation_RZSM
    ravel_init = 0
    ravel_len = int(ravel_length/sm_SubX_file.realization.shape[0])
    ravel_start = ravel_len
    for model in range(sm_SubX_file.realization.shape[0]):
        #ravel the first model, add to correlation_RZSM index column 1
        #index column 0 will be smerge
        
        
        correlation_RZSM[ravel_init:ravel_start,0] = sm_SMERGE_file.RZSM.values.ravel()
        correlation_RZSM[ravel_init:ravel_start,1] = sm_SubX_file.READ_DESCRIPTION_ATTRS[model,:,:,:].values.ravel()
        
        ravel_init += ravel_len
        ravel_start += ravel_len
    
    no_nan = correlation_RZSM[(~np.isnan(correlation_RZSM[:,0]) & (~np.isnan(correlation_RZSM[:,1])))]
    x,y = no_nan[:,0], no_nan[:,1]
    
    return(x,y,mean_x,mean_y)
    
#Plot heatmaps and scatterplots

try:
    plt.imread(f'{image_dir}/ETo_MSE_lagged_average.tif')
    plt.imread(f'{image_dir}/ETo_heatmap_lagged_average.tif')
    plt.imread(f'{image_dir}/ETo_MSE_lagged_average_mean.tif')
    plt.imread(f'{image_dir}/ETo_heatmap_lagged_average_mean.tif')
    
except FileNotFoundError:
    print('Finding MSE and heatmap for ETo lagged average ensemble.')    
    x,y,mean_x,mean_y = flatten_ETO()    
    
    #Plot ETo MSE
    max_ETo_val = 30
    
    x_line = np.linspace(0,max_ETo_val,100)
    y_line = np.linspace(0,max_ETo_val,100)
    
    try:
        plt.imread(f'{image_dir}/ETo_MSE_lagged_average.tif')
    except FileNotFoundError:
        ax = sns.scatterplot(x = x, y = y,s=0.25)
        ax.set_title("SubX ETo vs. SMERGE ETo")
        ax.set_xlabel("SMERGE (m3/m3)")
        ax.set_ylabel("SubX (m3/m3)")
        ax.set(xlim=(0,max_ETo_val))
        ax.set(ylim=(0,max_ETo_val))
        plt.plot(x_line,y_line,color='r')
        ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=9) #add text
        plt.savefig(f'{image_dir}/ETo_MSE_lagged_average.tif', dpi=300)
        plt.close()
    
    try:
        plt.imread(f'{image_dir}/ETo_heatmap_lagged_average.tif')
    except FileNotFoundError:
        #Plot heatmap
        heatmap, xedges, yedges = np.histogram2d(x,y, bins=100)
        extent = [0, yedges[-1], 0, yedges[-1]]
        plt.imshow(heatmap.T, extent=extent, origin='lower')
        plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(x,y),5)}", fontsize=14, color = 'red')
        plt.xlabel('SMERGE ETo (m3/m3)')
        plt.ylabel('SubX ETo (m3/m3)')
        plt.savefig(f'{image_dir}/ETo_heatmap_lagged_average.tif',dpi=300)
        plt.close()
        
        
    try:
        plt.imread(f'{image_dir}/ETo_MSE_lagged_average_mean.tif')
    except FileNotFoundError:
        ax = sns.scatterplot(x = mean_x, y = mean_y,s=0.25)
        ax.set_title("SubX ETo vs. SMERGE ETo")
        ax.set_xlabel("SMERGE (m3/m3)")
        ax.set_ylabel("SubX (m3/m3)")
        ax.set(xlim=(0,max_ETo_val))
        ax.set(ylim=(0,max_ETo_val))
        plt.plot(x_line,y_line,color='r')
        ax.text(0.01, 0.94,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=9) #add text
        plt.savefig(f'{image_dir}/ETo_MSE_lagged_average_mean.tif', dpi=300)
        plt.close()
        
         
    try:
        plt.imread(f'{image_dir}/ETo_heatmap_lagged_average_mean.tif')
    except FileNotFoundError:
    #Plot heatmap
        heatmap, xedges, yedges = np.histogram2d(mean_x, mean_y, bins=100)
        extent = [0, yedges[-1], 0, yedges[-1]]
        plt.imshow(heatmap.T, extent=extent, origin='lower')
        plt.text(0.10, 0.75,f"MSE: {np.round(mean_squared_error(mean_x,mean_y),5)}", fontsize=14, color = 'red')
        plt.xlabel('SMERGE ETo (m3/m3)')
        plt.ylabel('SubX ETo (m3/m3)')
        plt.savefig(f'{image_dir}/ETo_heatmap_lagged_average_mean.tif',dpi=300)
        plt.close()        
# var_OUT.READ_DESCRIPTION_ATTRS.values


#could get the mask to work here. May play around with it later.
# if var == 'msro':
#     #Save a dataset to have a soil moisture CONUS mask just like SMERGE
#     smerge_loc = f'{home}/Data/SMERGE_SM/Raw_data/'
#     SMERGE_SM = xr.open_dataset(f'{smerge_loc}/smerge_sm_merged_remap.nc4')
    
    
    
#     #Add new dimenesions and copy the first
#     SMERGE_SM_arr = np.array(SMERGE_SM.RZSM)
#     # SMERGE_SM_arr
#     # var_SM_mask[:,:,:,:] = SMERGE_SM_arr[:,:,:]
#     # var_SM_mask[var_SM_mask == 0.0] = 'nan' 
    
#     # var_where  = np.where(var_SM_mask[:,:,:,:]  != 0.0 )
    
#     # # var_SM_mask_where = np.where(var_SM_mask != 0.0)
#     # lagged_average_xr_mask = lagged_average_xr[var_where] 
#     # lagged_average_xr_mask1 = np.reshape(lagged_average_xr_mask, lagged_average_xr.shape)

### END mask test
'''Checking 2 different files to see if I'm actually getting the average

Uncomment for your own to self test and change the directory'''


# f1 = sorted(glob(f'{var}*.nc'))[0]
# f2 = sorted(glob(f'{var}*.nc'))[1]
# output_check = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/lagged_ensemble_by_variable/test1/mrso_lagged_average.nc4')

# #File1
# f1 = xr.open_dataset(f1)
# #File2
# f2 = xr.open_dataset(f2)
# #Check the date values
# f1.S.values[0] #1/11/1999 - contains values
# f1.S.values[1] #1/10/1999 - no values
# f2.S.values[1]

# #Get same days for each file... 5 day lag between these two files
# f1.mrso[0,0,5,0,0].values
# f2.mrso[0,0,0,0,0].values

# avg = (np.add(f1.mrso[0,0,5,0,0].values ,f2.mrso[0,0,0,0,0].values))/2.0 
# #convert to array 
# f1ar = np.array(f1.mrso)
# f2ar = np.array(f2.mrso)

# #Pick 6th day of the end time series.... Day 0 has no values
# print(f'Average test is complete and the results and the same? : {var_TEST[0,6,0,0] == avg}')

# #Correct calculation
