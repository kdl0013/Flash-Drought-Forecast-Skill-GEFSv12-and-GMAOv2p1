#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 09:06:40 2022

@author: kdl

Testing how to make a lagged ensemble average away from Cheyenne cluster


It appears that I'm making too large of a dimension, save each model into it's 
own file. Making a 1000 x 4 x 46 x 1045 x ? x ? is just too big
"""

import xarray as xr
import numpy as np
from time import sleep
import os
from glob import glob
from shutil import copyfile, copy
import pandas as pd
import datetime
from collections import OrderedDict
from multiprocessing import Pool, Manager, Process, Array
from time import time, sleep
from numba import njit, prange, float64, int64, vectorize

model = 'GMAO'
# variables = ['huss', 'dswrf','mrso','tas', 'uas', 'vas','tdps','pr','cape']


home_dir = f'/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/{model}/subset'
# num_processors = 5

os.chdir(home_dir)
os.getcwd()
var = 'mrso'

#get info for future dimension creation later
def inspect_files(var, date_start_val=-13, date_end_val=-3):
    '''Pre-processing, open the first file and inspect contents, save date list
    Inputs: 
        var = variable, date_start_val is the index location of date in file name
    Outputs: 
          list of dates in directory (date_list) and a single file to inspect (varFile) 
          and a count of files for that variable in directory. 
    '''
    count=0
    date_list = []
    for file in sorted(glob(f'{var}*.nc')):
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


date_list, varFile, count_files, daily_dates, varFileTEST = inspect_files(var=var)
#Convert date_list to a daily time series



varFile.info()

#Get info about date list and the actual number of days
number_days = pd.to_datetime(date_list[-1]).toordinal() - pd.to_datetime(date_list[0]).toordinal()
print('######### Number of Days in Final Time Series (Adding 44 days to end file)')
print(f'Total number of days in entire time series (plus lead days) is {number_days+varFile.L.shape[0]}')
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
dim_0 = 0
dim_1 = 1
for file in sorted(glob(f'{var}*.nc')[0:9]):
    print(f'File date says:  {file[-13:-3]}')
    file_o = xr.open_dataset(file)
    print(f'Number of missing values : {str(np.count_nonzero(np.isnan(file_o.mrso[dim_0,:,:,:,:].values)))}' + f'   Values are in S[{dim_0}] dim   Day {str(file_o.S.values[dim_0])[0:10]}') 
    print(f'Number of missing values : {str(np.count_nonzero(np.isnan(file_o.mrso[dim_1,:,:,:,:].values)))}' + f'    Values are in S[{dim_1}] dim   Day {str(file_o.S.values[dim_1])[0:10]}') 

#You can visually see that the first dimension has no np.nan, but S[1] dimension does
#%%Step 2, open a file and get the dates for each S dimension as well as create
#more dates based on lead time. 
#TODO: delete next line of code after testing
varFile=varFileTEST
#Stack a single file grid and make a copy to create a new empty dataset
varFile = varFile.stack(grid = ['Y','X'])

#TODO for each, model realization and each lead time run a loop
varFile = (varFile.mrso[0,0,0,:]).to_dataset()
varFile = varFile.drop(labels = ['L','M','S'])
#Add 2nd dimension for number of files

#File size list for all files, keeps up with index of file
varFile = expand_dimensions(varFile, file_size_list = np.array(range(0,count_files)))

#Add 3rd dimension for time series dates
varFile = expand_dimensions(varFile, time_series_dates = pd.to_datetime(daily_dates))

# #Add 4th dimension for model realization (don't do this, too much memory)
# varFile = expand_dimensions(varFile, realization = range(0,(file_o.M.shape[0])))

#Add 5th dimension for lead time, must +1 to account for the files initialization
#day that doesnt' contain data. lead 0 is actually the next days forecast
#Need this lead time dimension to stack properly
varFile = expand_dimensions(varFile, lead = range(0,(file_o.L.shape[0])+1))


varFile.mrso[:,:,:,:] = (0)
varFile= varFile.astype('float64')
varFile = varFile.unstack('grid')

var_TEST = varFile.copy()
#Output shape, use var_TEST as the output for the function
var_TEST.mrso.shape
#convert to numpy array
var_TEST = np.array(var_TEST.to_array())
var_TEST = var_TEST.squeeze()

#%%

'''I believe this function is working properly. 

Input: 
    File names
Output:
    var_TEST array which contains all the values of for each stacked dimesion
    from all SubX files. This will help to construct a soil moisture percentile
    climatology to include all values of each model and day  
'''

for file in sorted(glob('*.nc')):
    
    try:
        #Convert to numpy array, but keep date information
        def file_setup(file):
            #Read file
            open_f = xr.open_dataset(file)
            #Get date information
            beginning_start_date_ordinal = pd.to_datetime(list(sorted(glob('*.nc')))[0][-13:-3], format= '%Y/%m/%d').toordinal()
            beginning_start_date = pd.to_datetime(list(sorted(glob('*.nc')))[0][-13:-3])
            beginining_index = date_list.index((list(sorted(glob('*.nc')))[0][-13:-3]))
            
            #To keep up with file number within dimensions
            file_num = list(sorted(glob('*.nc'))).index(file) #Number of files. Each file in its own dimension (file_size_list)
            
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
            open_f_arr = open_f_arr[0,:,:,:,:]
            #Reshape
            # open_f_arr = np.reshape(open_f_arr,(open_f_arr.shape[0],open_f_arr.shape[2],open_f_arr.shape[1]))
            #Now we can maybe add Numba to fill in the missing values
            vector = np.vectorize(np.float)
            open_f_arr = vector(open_f_arr)
            
            return(open_f_arr, int(open_f_lead_0_ordinal), \
                   int(beginning_start_date_ordinal), int(file_num))
        
        open_f_arr, open_f_lead_0_ordinal, beginning_start_date_ordinal, file_num = file_setup(file)
    
        
        @vectorize([int64(int64,int64,int64,int64)])
        def lagged_ensemble_preprocess(open_f_lead_0_ordinal, beginning_start_date_ordinal, file_num, model):
            global open_f_arr
            global var_TEST
            # for model in prange(open_f_arr.shape[0]):
                # print(model)
            for lead in range(open_f_arr.shape[1]):
                lead_final_ordinal = open_f_lead_0_ordinal + lead
            # date_count = 0 #For understanding the starting date of the file
                for Yx in range(open_f_arr.shape[2]):   
                    # print(grid_stack)
                    for Xx in range(open_f_arr.shape[3]):
                
                    # print(open_f.dswrf[0,model,lead,grid_stack].data)
                    #To get the time_series_dates dimension data ordered correctly:
                    #var_TEST 3rd dimension data output will contain the index value of the
                    #lead_0 ordinal date when compared to the first file of the entire series' date
                    
                    #var_TEST output 
                        var_TEST[file_num,(lead_final_ordinal - beginning_start_date_ordinal),lead, Yx, Xx] = open_f_arr[model,lead, Yx, Xx] #append each lead to file_size_list dim
    
        model = 0
        #Run function
        lagged_ensemble_preprocess(open_f_lead_0_ordinal, beginning_start_date_ordinal, file_num, model)

    except TypeError:
        continue

#%%
#Convert to an xarray object
var_OUT = xr.Dataset(
    data_vars = dict(
        mrso = (['file_size_list','time_series_dates','lead','Y','X'], var_TEST),
    ),
    coords = dict(
        X = varFile.X.values,
        Y = varFile.Y.values,
        file_size_list = varFile.file_size_list.values,
        time_series_dates = daily_dates,
        lead = varFile.lead.values,
    ),
    attrs = dict(
        Description = 'Vertically integrated Soil Moisture'),
)

#Save as a netcdf for later processing
var_OUT.to_netcdf(path = 'Soil_moisture_full_dataset_model_0.nc4', mode ='w')


#%%




'''OLD CODE, multiprocessing that didn't work because the global interpreter lock
won't allow me to edit the same file while another file is working on it.'''




#This doesn't work, needed to vectorize
# @njit
# def lagged_ensemble_preprocess(open_f_arr,open_f_lead_0_ordinal, beginning_start_date_ordinal, file_num):
#     for model in prange(open_f_arr.shape[0]):
#         for grid_stack in range(open_f_arr.shape[1]):
#         # date_count = 0 #For understanding the starting date of the file
#             for lead in range(open_f_arr.shape[2]):    
#                 lead_final_ordinal = open_f_lead_0_ordinal + lead
            
#                 # print(open_f.dswrf[0,model,lead,grid_stack].data)
#                 #To get the time_series_dates dimension data ordered correctly:
#                 #var_TEST 3rd dimension data output will contain the index value of the
#                 #lead_0 ordinal date when compared to the first file of the entire series' date
                
#                 #var_TEST output 
#                 var_TEST[file_num,grid_stack,model,(lead_final_ordinal - beginning_start_date_ordinal),lead] = 1 #append each lead to file_size_list dim






# open_f_arr.shape
# var_TEST.shape

# '''Currently this code block will work, but it sucks and its slow, 
# testing Numba now... and maybe dask.... I'm going with Numba'''

# #I could also convert the function to a numba function... just need to re-order the 
# #variables... it would be probably faster

# #TODO: Figure out why multiprocessing won't input data into the file....
# #Will have to switch to a loop which would take forever...
# #Other option is to convert it all into a numba process and then convert back
# #to an xarray object... that may be the simplest.

# out_array = np.empty(var_TEST.mrso.shape)

# def lagged_ensemble_preprocess(file):
#     '''Combine all files into a single file as a lagged ensemble, but do not take
#     the average. This will be used for soil moisture files to get a soil moisture
#     percentile that captures all of the models values to increase the distribution
#     total values to draw from.
    
#     Inputs: List of files (entered as self in multiprocessin function below)
#     Output: Single file with dimensions:
#         (num_files) x 1593 (XY-grid-stacked) x 4 (realizations) x 9 (date_list) x daily time series (as.Date)
#     '''
#     print(f'Starting {file}.  ')
    
#     global var_TEST
# #Loop works, but it would be very very slow (about )    
# # for file in list(sorted(glob(f'{var}*.nc'))):
#     beginning_start_date_ordinal = pd.to_datetime(list(sorted(glob('*.nc')))[0][-13:-3], format= '%Y/%m/%d').toordinal()
#     beginning_start_date = pd.to_datetime(list(sorted(glob('*.nc')))[0][-13:-3])
#     beginining_index = date_list.index((list(sorted(glob('*.nc')))[0][-13:-3]))
    
#     #To keep up with file number within dimensions
#     file_num = list(sorted(glob('*.nc'))).index(file) #Number of files. Each file in its own dimension (file_size_list)
    
#     #All opened files are files that have been initialized, pre-processing has
#     #already occurred to remove non-initialized days
#     open_f = xr.open_dataset(file).stack(grid = ['Y','X'])
    
#     #Get date info about opened, intialized file
#     open_f_lead_0_date = str(pd.to_datetime(open_f.S.values[0], format= '%Y/%m/%d'))[0:10] #Character string of lead 0 day (day after initialization)
#     open_f_lead_0_minus_1_date = str(pd.to_datetime(open_f.S.values[1], format= '%Y/%m/%d'))[0:10]
    
#     open_f_lead_0_ordinal = pd.to_datetime(open_f_lead_0_date, format= '%Y/%m/%d').toordinal()
#     open_f_lead_0_minus_1_ordinal = pd.to_datetime(open_f_lead_0_minus_1_date, format= '%Y/%m/%d').toordinal()
    
#     open_f_lead_0_minus_1_index = file_num  #Same as file n
    
#     for grid_stack in range(open_f.grid.shape[0]):
#         for model in range(open_f.M.shape[0]):
#         # date_count = 0 #For understanding the starting date of the file
#             for lead in range(open_f.L.shape[0]):    
#                 lead_final_ordinal = open_f_lead_0_ordinal + lead
            
#                 # print(open_f.dswrf[0,model,lead,grid_stack].data)
#                 #To get the time_series_dates dimension data ordered correctly:
#                 #var_TEST 3rd dimension data output will contain the index value of the
#                 #lead_0 ordinal date when compared to the first file of the entire series' date
                
#                 #var_TEST output 
#                 var_TEST.mrso[file_num,grid_stack,model,(lead_final_ordinal - beginning_start_date_ordinal),lead] = open_f.mrso[0,model,lead,grid_stack].values #append each lead to file_size_list dim

#     print(f'File {file_num} of {len(glob("*.nc"))} completed.')
#     # return(var_TEST)


# if __name__ == '__main__':
    
#     # arr = Array()
#     # p = Process(target = lagged_ensemble_preprocess, args = list(sorted(glob(f'{var}*.nc'))), )
    
    
#     p = Pool(5)
#     start_t = time()
    
#     # # global var_TEST
    
#     # # var_TEST = var_TEST
#     p.map(lagged_ensemble_preprocess, list(sorted(glob(f'{var}*.nc'))))
    
#     end_t = time()
#     len_files = len(sorted(glob(f'{var}*.nc')))
#     print(f'{(end_t - start_t) / 60} minutes for {len_files} files.')
#     var_TEST.mrso.values
#     np.nansum(open_f.mrso.values)
#     np.nansum(var_TEST.mrso.values)
#     var_TEST.mrso.shape
#     var_TEST.mrso[0,:,0,1,:].shape
#     var_TEST.mrso[0,:,0,1,:].values
#     var_TEST.mrso[0,:,:,:,:].values.sum()
#     var_TEST_save = var_TEST.unstack('grid')
#     var_TEST_save.to_netcdf(path = f'{home_dir}/test.nc')



