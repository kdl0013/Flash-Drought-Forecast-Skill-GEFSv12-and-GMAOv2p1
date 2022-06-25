#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 18:02:41 2022

@author: kdl
"""


import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
import refet

dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
gridMET_dir = f'{dir1}/Data/gridMET'
fileOUT_gridMET = f'{gridMET_dir}/ETo_anomaly_gridMET_merged.nc'

SMERGE_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
fileOUT_SMERGE = f'{SMERGE_dir}/RZSM_anomaly_SMERGE_merged.nc'

#gridMET ETo
open_f = xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
anomaly_f = xr.zeros_like(open_f)
anomaly_date = pd.DataFrame(open_f.ETo_gridmet.day)
anomaly_date.index = pd.to_datetime(anomaly_date.iloc[:,0])
anomaly_date_list = list(anomaly_date.iloc[:,0])

#SMERGE RZSM
open_rzsm = xr.open_dataset(f'{SMERGE_dir}/smerge_sm_merged_remap.nc4')
anomaly_r = xr.zeros_like(open_rzsm)
anomaly_date_r = pd.DataFrame(open_rzsm.RZSM.time)
anomaly_date_r.index = pd.to_datetime(anomaly_date_r.iloc[:,0])
anomaly_date_r_list = list(anomaly_date_r.iloc[:,0])

#gridmet anomalies
try:
    xr.open_dataset(fileOUT_gridMET) #see if file is already created
    print(f'Already created gridMET ETo anomaly file into {gridMET_dir}.')
except FileNotFoundError:
    
    
    for lat in range(open_f.ETo_gridmet.shape[1]):
        print(f'Working on lat {lat} out of {np.max(range(open_f.ETo_gridmet.shape[1]))} for ETo.')
        for lon in range(open_f.ETo_gridmet.shape[2]):
    
            # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
            #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
            #Only need to look at 1 years worth of data because we are adding 
            for day in range(7,377):
                
                #Don't re-run additionaly data if its already in the final file
                if anomaly_f.ETo_gridmet.isel(lat=lat,lon=lon, day=day).values != 0:
                    pass
                else:
                    #find the date of the current step
                    count=0
                    new_day_val = open_f.day[day-7].values 
                    mean_dict = {}
                    #Total of 17 years worth of data from SubX
                    while count!=17:
                        count+=1
                        weekly_mean = np.nanmean(open_f.ETo_gridmet.isel(lat=lat,lon=lon).sel(day=slice(new_day_val-np.timedelta64(7,'D'),new_day_val)))
                        #add to a list
                        mean_dict[f'{new_day_val}']=weekly_mean
                        #if leap year, add 1 year to loop through all years
                        if pd.to_datetime(new_day_val).year % 4 ==0:
                            new_day_val = new_day_val+np.timedelta64(366,'D')
                        else:
                            new_day_val = new_day_val+np.timedelta64(365,'D')
                    
                    #Now that we have the mean value for each year, find the mean of the means and subtract to make anomaly
                    mean_vals = np.nanmean([i for i in mean_dict.values()])
                    #Now subtract the mean from all values
                    for i in mean_dict.keys():
                        mean_dict[i] = mean_dict[i] - mean_vals
                        
                    #Now add to the empty anomly_f file
                    for i in mean_dict.keys():
                        _date = pd.to_datetime(i)
                        #find the index in anomaly_date_list
                        index_val = anomaly_date_list.index(_date)
                        #add to file
                        anomaly_f.ETo_gridmet[index_val, lat,lon] = mean_dict[i]
                        
            # '''Don't know why, I have to re-run this again to make all values appear'''
            # for lat in range(open_f.ETo_gridmet.shape[1]):
            #     for lon in range(open_f.ETo_gridmet.shape[2]):
                 
            #         for day in range(360,374):
                        
            #             #Don't re-run additionaly data if its already in the final file
            #             if anomaly_f.ETo_gridmet.isel(lat=lat,lon=lon, day=day).values != 0:
            #                 pass
            #             else:
            #                 #find the date of the current step
            #                 count=0
            #                 new_day_val = open_f.day[day-7].values 
            #                 mean_dict = {}
            #                 #Total of 17 years worth of data from SubX
            #                 while count!=17:
            #                     count+=1
            #                     weekly_mean = np.nanmean(open_f.ETo_gridmet.isel(lat=lat,lon=lon).sel(day=slice(new_day_val-np.timedelta64(7,'D'),new_day_val)))
            #                     #add to a list
            #                     mean_dict[f'{new_day_val}']=weekly_mean
            #                     #if leap year, add 1 year to loop through all years
            #                     if pd.to_datetime(new_day_val).year % 4 ==0:
            #                         new_day_val = new_day_val+np.timedelta64(366,'D')
            #                     else:
            #                         new_day_val = new_day_val+np.timedelta64(365,'D')
                            
            #                 #Now that we have the mean value for each year, find the mean of the means and subtract to make anomaly
            #                 mean_vals = np.nanmean([i for i in mean_dict.values()])
            #                 #Now subtract the mean from all values
            #                 for i in mean_dict.keys():
            #                     mean_dict[i] = mean_dict[i] - mean_vals
                                
            #                 #Now add to the empty anomly_f file
            #                 for i in mean_dict.keys():
            #                     _date = pd.to_datetime(i)
            #                     #find the index in anomaly_date_list
            #                     index_val = anomaly_date_list.index(_date)
            #                     #add to file
            #                     anomaly_f.ETo_gridmet[index_val, lat,lon] = mean_dict[i]
            
        anomaly_f.to_netcdf(path = fileOUT_gridMET, mode ='w', engine='scipy')
        anomaly_f.close()

#%%
'''SMERGE anomaly (use same methodology as gridMET ETo anomaly). This is the same
methodology as SubX as well for anomaly calcualation. This also keeps the distribution
within a select number of years'''

#smerge anomalies
try:
    xr.open_dataset(fileOUT_SMERGE) #see if file is already created
    print(f'Already created SMERGE RZSM anomaly file into {SMERGE_dir}.')

except FileNotFoundError:
    # lat=10
    # lon=10
    
    for lat in range(open_rzsm.RZSM.shape[1]):
        print(f'Working on lat {lat} out of {np.max(range(open_rzsm.RZSM.shape[1]))} for RZSM.')
        for lon in range(open_rzsm.RZSM.shape[2]):
    
            # all_values = open_rzsm.RZSM.isel(Y=lat,X=lon)
            #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
            #Only need to look at 1 years worth of data because we are adding 
            for day in range(7,377):
                
                #Don't re-run additionaly data if its already in the final file
                if anomaly_r.CCI_ano.isel(Y=lat,X=lon, time=day).values != 0:
                    pass
                else:
                    #find the date of the current step
                    count=0
                    new_day_val = open_f.day[day-7].values 
                    mean_dict = {}
                    #Total of 17 years worth of data from SubX
                    while count != 17:
                        count+=1
                        weekly_mean = np.nanmean(open_rzsm.RZSM.isel(X=lon,Y=lat).sel(time=slice(new_day_val-np.timedelta64(7,'D'),new_day_val)))
                        #add to a list
                        mean_dict[f'{new_day_val}']=weekly_mean
                        #if leap year, add 1 year to loop through all years
                        if pd.to_datetime(new_day_val).year % 4 ==0:
                            new_day_val = new_day_val+np.timedelta64(366,'D')
                        else:
                            new_day_val = new_day_val+np.timedelta64(365,'D')
                    
                    #Now that we have the mean value for each year, find the mean of the means and subtract to make anomaly
                    mean_vals = np.nanmean([i for i in mean_dict.values()])
                    #Now subtract the mean from all values
                    for i in mean_dict.keys():
                        mean_dict[i] = mean_dict[i] - mean_vals
                        
                    #Now add to the empty anomly_f file
                    for i in mean_dict.keys():
                        _date = pd.to_datetime(i)
                        #find the index in anomaly_date_list
                        index_val = anomaly_date_r_list.index(_date)
                        #add to file
                        anomaly_r.CCI_ano[index_val, lat,lon] = mean_dict[i]
                        
                      
        # for lat in range(open_rzsm.RZSM.shape[1]):
        #     for lon in range(open_rzsm.RZSM.shape[2]):
             
        #         for day in range(360,374):
                    
        #             #Don't re-run additionaly data if its already in the final file
        #             if anomaly_r.CCI_ano.isel(Y=lat,X=lon, time=day).values != 0:
        #                 pass
        #             else:
        #                 #find the date of the current step
        #                 count=0
        #                 new_day_val = open_f.day[day-7].values 
        #                 mean_dict = {}
        #                 #Total of 17 years worth of data from SubX
        #                 while count!=17:
        #                     count+=1
        #                     weekly_mean = np.nanmean(open_rzsm.RZSM.isel(X=lon,Y=lat).sel(time=slice(new_day_val-np.timedelta64(7,'D'),new_day_val)))
        #                     #add to a list
        #                     mean_dict[f'{new_day_val}']=weekly_mean
        #                     #if leap year, add 1 year to loop through all years
        #                     if pd.to_datetime(new_day_val).year % 4 ==0:
        #                         new_day_val = new_day_val+np.timedelta64(366,'D')
        #                     else:
        #                         new_day_val = new_day_val+np.timedelta64(365,'D')
                        
        #                 #Now that we have the mean value for each year, find the mean of the means and subtract to make anomaly
        #                 mean_vals = np.nanmean([i for i in mean_dict.values()])
        #                 #Now subtract the mean from all values
        #                 for i in mean_dict.keys():
        #                     mean_dict[i] = mean_dict[i] - mean_vals
                            
        #                 #Now add to the empty anomly_f file
        #                 for i in mean_dict.keys():
        #                     _date = pd.to_datetime(i)
        #                     #find the index in anomaly_date_list
        #                     index_val = anomaly_date_list.index(_date)
        #                     #add to file
        #                 anomaly_r.CCI_ano[index_val, lat,lon] = mean_dict[i]
                        
    anomaly_r.to_netcdf(path = fileOUT_SMERGE, mode ='w', engine='scipy')
    anomaly_r.close()