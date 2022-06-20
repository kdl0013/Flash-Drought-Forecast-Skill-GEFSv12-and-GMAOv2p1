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
from multiprocessing import Pool


dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
gridMET_dir = f'{dir1}/Data/gridMET'
fileOUT = f'{gridMET_dir}/ETo_anomaly_gridMET_merged.nc'

open_f = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/gridMET/ETo_gridMET_merged.nc')
anomaly_f = xr.zeros_like(open_f)
anomaly_date = pd.DataFrame(open_f.ETo_gridmet.day)
anomaly_date.index = pd.to_datetime(anomaly_date.iloc[:,0])
anomaly_date_list = list(anomaly_date.iloc[:,0])



for lat in range(open_f.ETo_gridmet.shape[1]):
    print(f'Working on lat {lat}.')
    for lon in range(open_f.ETo_gridmet.shape[2]):

        all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
        #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
        #Only need to look at 1 years worth of data because we are adding 
        for day in range(7,372):
            
            #Don't re-run additionaly data if its already in the final file
            if anomaly_f.ETo_gridmet.isel(lat=lat,lon=lon, day=day).values != 0:
                pass
            else:
                #find the date of the current step
                count=0
                new_day_val = open_f.day[day-7].values 
                mean_dict = {}
                #Total of 17 years worth of data from 
                while count!=17:
                    count+=1
                    weekly_mean = open_f.ETo_gridmet.isel(lat=lat,lon=lon).sel(day=slice(new_day_val-np.timedelta64(7,'D'),new_day_val)).mean()
                    #add to a list
                    mean_dict[f'{new_day_val}']=weekly_mean
                    #if leap year
                    if pd.to_datetime(new_day_val).year % 4 ==0:
                        new_day_val = new_day_val+np.timedelta64(366,'D')
                    else:
                        new_day_val = new_day_val+np.timedelta64(365,'D')
                
                #Now that we have the mean value for each year, find the mean of the means and subtract to make anomaly
                mean_vals = np.mean([i for i in mean_dict.values()])
                #Now subtract the mean from all values
                for i in mean_dict.keys():
                    mean_dict[i] = mean_dict[i] - mean_vals
                    
                #Now add to the empty anomly_f file
                for i in mean_dict.keys():
                    _date = pd.to_datetime(i)
                    #find the index in anomaly_date_list
                    index_val = anomaly_date_list.index(_date)
                    #add to file
                    anomaly_f.ETo_gridmet[index_val, lat,lon] = mean_dict[i].values

anomaly_f.to_netcdf(path = fileOUT, mode ='w', engine='scipy')
anomaly_f.close()