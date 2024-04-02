#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 18:02:41 2022

@author: kdl
"""

import xarray as xr
import numpy as np
import pandas as pd

# dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
dir1 = 'main_dir'
gridMET_dir = f'{dir1}/Data/gridMET'
fileOUT_gridMET = f'{gridMET_dir}/ETo_anomaly_gridMET_merged.nc'
fileOUT_gridMET_mean =f'{gridMET_dir}/ETo_anomaly_gridMET_mean.nc'

SMERGE_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
fileOUT_SMERGE = f'{SMERGE_dir}/RZSM_anomaly_SMERGE_merged.nc'
fileOUT_SMERGE_mean = f'{SMERGE_dir}/RZSM_anomaly_SMERGE_mean.nc'

global anomaly_range
anomaly_range=42 #choose 42 days on either side of date to create anomaly

#gridMET ETo
open_f = xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
anomaly_f = xr.zeros_like(open_f)
anomaly_f_mean = xr.zeros_like(open_f)

anomaly_date = pd.DataFrame(open_f.ETo_gridmet.day)
anomaly_date.index = pd.to_datetime(anomaly_date.iloc[:,0])
anomaly_date_list = list(anomaly_date.iloc[:,0])

# #get date of list
# day_list=pd.to_datetime(anomaly_f.day.values)
# out_day_list = []
# for day_ in day_list:
#     period = pd.Period(day_, freq='D')
#     doy = period.day_of_year
#     out_day_list.append(doy)

# anomaly_f = anomaly_f.assign_coords(day=out_day_list)

# #Perform a 7-day rolling mean
# open_f = open_f.assign_coords(day=out_day_list)
# open_f.ETo_gridmet[:,10,10].values

#SMERGE RZSM
open_rzsm = xr.open_dataset(f'{SMERGE_dir}/smerge_sm_merged_remap.nc4')
anomaly_r = xr.zeros_like(open_rzsm)
anomaly_r_mean = xr.zeros_like(anomaly_r)
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

                #Calculate anomaly based on a 42 day window on both sides of the 
                #day (based on Julian day only)
                #find the date of the current step
                day_val=open_f.day[day].to_numpy()
          
                yearly_average = {}
                new_day_val = day_val
                all_values = []

                #Total of 17 years worth of data from SubX
                for year_ in np.arange(1999,2017):

                    weekly_mean = np.nanmean(open_f.ETo_gridmet.isel(lat=lat,lon=lon).sel(day=slice(day_val-np.timedelta64(7,'D'),day_val)))
                    #add to a list
                    # all_values.append(weekly_mean)
                    yearly_average[f'{day_val}']=weekly_mean
                    #if leap year, add 1 year to loop through all years
                    if pd.to_datetime(day_val).year % 4 ==0:
                        day_val = day_val+np.timedelta64(366,'D')
                    else:
                        day_val = day_val+np.timedelta64(365,'D')
                
                yearly_avg_list = [i for i in yearly_average.values()]
                #Now choose days 42 days on either side of each year and append to a single file
                day_val=open_f.day[day].to_numpy()
                for year_ in np.arange(1999,2017):
                    if pd.to_datetime(day_val).year % 4 ==0:
                        day_val = day_val+np.timedelta64(366,'D')
                    else:
                        day_val = day_val+np.timedelta64(365,'D')
                    all_values.append(list(open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()))
                    # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                
                all_values.append(yearly_avg_list)
                #flatten the list of lists, find mean
                mean_ = np.nanmean([x for xs in all_values for x in xs])
                
                #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
                
                
                #subtract mean from each file in yearly_average dictionary (anomaly)
                for date_ in yearly_average.keys():
                    yearly_average[date_] = yearly_average[date_] - mean_
   
                #Now add to the empty anomly_f file
                for i in yearly_average.keys():
                    _date = pd.to_datetime(i)
                    #find the index in anomaly_date_list
                    try:
                        index_val = anomaly_date_list.index(_date)
                        anomaly_f.ETo_gridmet[index_val, lat,lon] = yearly_average[i]
                    except ValueError:
                        pass
                    
                    try:
                        index_val = anomaly_date_list.index(_date)
                        anomaly_f_mean.ETo_gridmet[index_val, lat,lon] = mean_
                    except ValueError:
                        pass

            
        anomaly_f.to_netcdf(path = fileOUT_gridMET, mode ='w', engine='scipy')
        anomaly_f_mean.to_netcdf(path = fileOUT_gridMET_mean, mode ='w', engine='scipy')

        anomaly_f.close()

#%%
'''SMERGE anomaly (use same methodology as gridMET ETo anomaly). This is the same
methodology as SubX as well for anomaly calcualation. This also keeps the distribution
within a select number of years'''


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
                
                day_val=open_f.day[day].to_numpy()
          
                yearly_average = {}
                new_day_val = day_val
                all_values = []

                #Total of 17 years worth of data from SubX
                for year_ in np.arange(1999,2017):

                    weekly_mean = np.nanmean(open_rzsm.RZSM.isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
                    #add to a list
                    # all_values.append(weekly_mean)
                    yearly_average[f'{day_val}']=weekly_mean
                    #if leap year, add 1 year to loop through all years
                    if pd.to_datetime(day_val).year % 4 ==0:
                        day_val = day_val+np.timedelta64(366,'D')
                    else:
                        day_val = day_val+np.timedelta64(365,'D')
                
                yearly_avg_list = [i for i in yearly_average.values()]
                #Now choose days 42 days on either side of each year and append to a single file
                day_val=open_f.day[day].to_numpy()
                for year_ in np.arange(1999,2017):
                    if pd.to_datetime(day_val).year % 4 ==0:
                        day_val = day_val+np.timedelta64(366,'D')
                    else:
                        day_val = day_val+np.timedelta64(365,'D')
                    all_values.append(list(open_rzsm.RZSM.sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))
                    # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
               
                all_values.append(yearly_avg_list)
                #flatten the list of lists, find mean
                mean_ = np.nanmean([x for xs in all_values for x in xs])
               
                #subtract mean from each file in yearly_average dictionary (anomaly)
                for date_ in yearly_average.keys():
                    yearly_average[date_] = yearly_average[date_] - mean_
   
                #Now add to the empty anomly_f file
                for i in yearly_average.keys():
                    _date = pd.to_datetime(i)
                    #find the index in anomaly_date_list
                    try:
                        index_val = anomaly_date_r_list.index(_date)
                        anomaly_r.CCI_ano[index_val, lat,lon] = yearly_average[i]
                    except ValueError:
                        pass
                    
                    try:
                        index_val = anomaly_date_list.index(_date)
                        anomaly_r_mean.CCI_ano[index_val, lat,lon] = mean_
                    except ValueError:
                        pass
                    
        anomaly_r.to_netcdf(path = fileOUT_SMERGE, mode ='w', engine='scipy')
        anomaly_r_mean.to_netcdf(path = fileOUT_SMERGE_mean, mode ='w', engine='scipy')

        anomaly_r.CCI_ano[index_val, lat,lon] = mean_
        anomaly_r.close()