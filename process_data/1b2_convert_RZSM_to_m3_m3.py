#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert GMAO soil moisture from kg/m2 to m3/m3

@author: kdl
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
from datetime import timedelta
import pandas as pd
from glob import glob
from multiprocessing import Pool


dir1 = 'main_dir'
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

# model='GMAO'
subX_dir = f'{dir1}/Data/SubX/{model}'
subX_SM_out_dir = f'{subX_dir}/SM_converted_to_m3_m3'

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc')):
        date_list.append(file[-13:-3])
    return(date_list)
        
init_date_list = return_date_list(var = 'mrso')    
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
        
        SubX_file = xr.open_dataset(f'{subX_dir}/mrso_GMAO_{_date}.nc4')
        var='mrso'
        
        def date_file_info(SubX_file):

            a_date_in= SubX_file.L.values
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

            return(a_date_out)
        
        julian_list = date_file_info(SubX_file)
        
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
            SM_SubX_m3_m3_value = (['S','model','lead','Y','X'], out_sm_SubX.mrso.values),
        ),
        coords = dict(
        S = out_sm_SubX.S.values,
        X = out_sm_SubX.X.values,
        Y = out_sm_SubX.Y.values,
        lead = julian_list,
        model = out_sm_SubX.M.values,
        ),
        attrs = dict(
        Description = 'SubX SM (kg/m2) converted to (m3/m3) based on profile depth.'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{subX_SM_out_dir}/SM_SubX_m3_m3_{_date}.nc', mode ='w')
        print(f'Completed SM_m3_m3_{_date}.')


'''I Moved this script, fix later to run... multiprocessing funcitons'''
def convert_SubX_to_m3_m3(_date,var):
