#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create evaporative demand.

Create 42 day spread for anomalies for MERRA soil moisture and evaporative demand.


@author: kdl
"""

import xarray as xr
import numpy as np
import pandas as pd
from metpy.units import units
import metpy.calc
import os

dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
MERRA_dir = f'{dir1}/Data/MERRA2'

evap='Priestley'
eto_merra_save_name = f'{MERRA_dir}/ETo_{evap}_merged.nc'

fileOUT_MERRA = f'{MERRA_dir}/ETo_anomaly_{evap}_MERRA_mean.nc'
fileOUT_MERRA_eto_anomaly = f'{MERRA_dir}/ETo_anomaly_{evap}_MERRA.nc'

fileOUT_MERRA_root = f'{MERRA_dir}/RZSM_anomaly_MERRA_mean.nc'
fileOUT_MERRA_root_anomaly = f'{MERRA_dir}/RZSM_anomaly_MERRA.nc'

fileOUT_MERRA_windspeed = f'{MERRA_dir}/windspeed_anomaly_MERRA_mean.nc'
fileOUT_MERRA_windspeed_anomaly = f'{MERRA_dir}/windspeed_anomaly_MERRA.nc'

fileOUT_MERRA_radiation = f'{MERRA_dir}/radiation_anomaly_MERRA_mean.nc'
fileOUT_MERRA_radiation_anomaly = f'{MERRA_dir}/radiation_anomaly_MERRA.nc'

fileOUT_MERRA_tasmax = f'{MERRA_dir}/tasmax_anomaly_MERRA_mean.nc'
fileOUT_MERRA_tasmax_anomaly = f'{MERRA_dir}/tasmax_anomaly_MERRA.nc'

fileOUT_MERRA_tasmin = f'{MERRA_dir}/tasmin_anomaly_MERRA_mean.nc'
fileOUT_MERRA_tasmin_anomaly = f'{MERRA_dir}/tasmin_anomaly_MERRA.nc'

fileOUT_MERRA_tas = f'{MERRA_dir}/tas_anomaly_MERRA_mean.nc'
fileOUT_MERRA_tas_anomaly = f'{MERRA_dir}/tas_anomaly_MERRA.nc'

fileOUT_MERRA_actual_vapor_pressure = f'{MERRA_dir}/actual_vapor_pressure_anomaly_MERRA_mean.nc'
fileOUT_MERRA_actual_vapor_pressure_anomaly = f'{MERRA_dir}/actual_vapor_pressure_anomaly_MERRA.nc'

CONUS_mask = xr.open_dataset(f'{dir1}/Scripts/CONUS_mask/NCA-LDAS_masks_SubX.nc4')

global anomaly_range
anomaly_range=42 #choose 42 days on either side of date to create anomaly

# #gridMET ETo
# open_f = xr.open_dataset(f'{MERRA_dir}/ETo_merged.nc4')
# anomaly_f = xr.zeros_like(open_f)
# anomaly_f_mean = xr.zeros_like(open_f)

# anomaly_date = pd.DataFrame(open_f.ETo_gridmet.day)
# anomaly_date.index = pd.to_datetime(anomaly_date.iloc[:,0])
# anomaly_date_list = list(anomaly_date.iloc[:,0])

open_rzsm = xr.open_dataset(f'{MERRA_dir}/RZSM_1m_column.nc4') #1 meter column
radiation = xr.open_dataset(f'{MERRA_dir}/radiation.nc4')
temp = np.subtract(xr.open_dataset(f'{MERRA_dir}/temperature.nc4'),273.15) #open, convert to Celsius
ptC = xr.open_dataset(f'{dir1}/Data/Priestley_Taylor_Makkink_evap_coeff_maps/pt_coeff_FINAL.nc') #priestley-taylor coefficient
elevation = xr.open_dataset(f'{dir1}/Data/elevation/elev_regrid.nc',decode_times=False)
vapor_pressure = xr.open_dataset(f'{MERRA_dir}/actual_vapor_pressure.nc4')


anomaly_r = xr.zeros_like(open_rzsm) #OUTPUT
#Get datesopen_f
anomaly_date_r = pd.DataFrame(open_rzsm[list(open_rzsm.keys())[0]].time)
anomaly_date_r.index = pd.to_datetime(anomaly_date_r.iloc[:,0])
anomaly_date_r_list = list(anomaly_date_r.iloc[:,0])

wind_humidity = xr.open_dataset(f'{MERRA_dir}/wind_humidity.nc4')
wind = wind_humidity.SPEED.to_dataset().rename(SPEED='windspeed')

def name_1(dict_):
    return(list(dict_.keys())[0])

def name_2(dict_):
    return(list(dict_.keys())[1])


#%% RZSM
#RZSM anomalies
def run_MERRA_RZSM(open_rzsm):
    anomaly_r = xr.zeros_like(open_rzsm) #OUTPUT
    anomaly_r = anomaly_r.drop('time_bnds')
    anomaly_r = anomaly_r.rename(RZMC = 'RZSM_mean')
    try:
        check_mean_file=xr.open_dataset(fileOUT_MERRA_root) #see if file is already created
        print(f'Already created MERRA RZSM anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for RZSM anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_2(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_2(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = fileOUT_MERRA_root, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_2(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('RZSM_anom')
                anomaly_create.to_netcdf(fileOUT_MERRA_root_anomaly, mode ='w', engine='scipy')
                anomaly_create.close()
                    
    
    return(0)

#%% ETo Priestley
def make_MERRA_ETo():
    try:
        save_output=xr.open_dataset(eto_merra_save_name)
    except FileNotFoundError:
        total_radiation = (radiation.SWNETSRF + radiation.LWGNET).to_dataset(name='net_radiation')
        
        #convert from W/m2 to MJ/m2/d
        mult_ = 0.0864
        convert_rad = np.multiply(total_radiation,mult_)
        #save output
        save_output = xr.zeros_like(total_radiation).rename(net_radiation='ETo')
        
        #need to only loop through days because ptC doesn't have any dates in the file
        print('Starting Priestley Taylor calculation')
        for day in range((temp.time.shape[0])):
            if day%100==0:
                print(f'Completed number {day} out of {max(range(temp.time.shape[0]))}')
            
            def priestley_taylor(temp1,convert_rad1,ptC1,elevation1):
                atm_pressure = 101.3 * ((293 - 0.0065 * elevation1)/293)**5.26 # kPa --- LOOKS FINE AFTER INSPECTION
                specific_heat = 0.001013 #MJ/kg C
                mle_ratio_weight_vapor = 0.622
                latent_heat_of_vaporization = 2.501 - (0.002361 * temp1) #MJ/kg, needs temp in Celsius
                
                '''gamma practically same results as old_gamma variable which uses only 0.000665 as the constant'''
                gamma = np.divide(np.multiply(atm_pressure,specific_heat), (mle_ratio_weight_vapor*latent_heat_of_vaporization)) #kPA
                # old_gamma = np.multiply(0.000665, atm_pressure) # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
                
                delta = 4098 * (0.6108 * np.exp(17.27 * temp1 / (temp1  + 237.3))) / (temp1  + 237.3)**2 # Slope saturated vapor pressure #kPA
                soil_heat_flux = np.zeros_like(delta)
                PET = (ptC1*(delta/(delta+gamma))*(convert_rad1-soil_heat_flux))/latent_heat_of_vaporization
    
                return PET
    
            #OUtput is saved in meters
            save_output.ETo[day,:,:]=priestley_taylor(temp1=temp.T2MMEAN[day,:,:],\
                              convert_rad1=convert_rad[list(convert_rad.keys())[0]][day,:,:],\
                                  ptC1=ptC[list(ptC.keys())[0]][:,:],\
                                      elevation1=elevation[list(elevation.keys())[0]][0,:,:])        # def priestley_taylor_new_equation(temp1,convert_rad1,ptC1,elevation1):
                
            
           
            # temp1=temp[list(temp.keys())[0]][day,10,10].values
            # convert_rad1=convert_rad[list(convert_rad.keys())[0]][day,10,10].values
            # ptC1=ptC[list(ptC.keys())[0]][10,10].values
            # elevation1=elevation[list(elevation.keys())[0]][0,10,10].values    
        
            # #OUtput is saved in meters
            # temp1=temp[list(temp.keys())[0]][day,:,:]
            # convert_rad1=convert_rad[list(convert_rad.keys())[0]][day,:,:]
            # ptC1=ptC[list(ptC.keys())[0]][:,:]
            # elevation1=elevation[list(elevation.keys())[0]][0,:,:]       # def priestley_taylor_new_equation(temp1,convert_rad1,ptC1,elevation1):
    
         # save_output = np.divide(save_output,1000)#now in mm (from old equation)
        
        save_output.to_netcdf(f'{eto_merra_save_name}',mode='w')
        #Create a yearly sum file for inspection
        os.system(f'rm {MERRA_dir}/eto_{evap}_monthly_mean.nc && cdo monmean {eto_merra_save_name} {MERRA_dir}/eto_{evap}_monthly_mean.nc')
        os.system(f'rm {MERRA_dir}/eto_{evap}_mean_annual.nc && cdo timmean -yearsum {eto_merra_save_name} {MERRA_dir}/eto_{evap}_mean_annual.nc')
        os.system(f'rm {MERRA_dir}/eto_{evap}_year_average.nc && cdo -yearsum {eto_merra_save_name} {MERRA_dir}/eto_{evap}_year_average.nc')
        
        return(save_output)

#%% Priestley anomaly
#ETo anomalies
def run_MERRA_ETo_anomaly():
    open_eto = xr.open_dataset(f'{eto_merra_save_name}')
    #Create empty files rzsm
    anomaly_e = xr.zeros_like(open_eto)
    anomaly_e_mean = xr.zeros_like(anomaly_e)

    try:
        anomaly_e=xr.open_dataset(fileOUT_MERRA) #see if file is already created
        print(f'Already created MERRA anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
    
        for lat in range(open_eto.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for ETo anomaly MERRA2.')
            for lon in range(open_eto.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_eto.time[day].to_numpy()
                        
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_eto.ETo.isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_eto.ETo.sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
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
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_e.ETo[index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_e.ETo[index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                                
            anomaly_e = anomaly_e.rename(ETo = 'ETo_mean')
            anomaly_e.ETo_mean[0:7,:,:] = anomaly_e.ETo_mean[0+366:7+366].values #account for the first year because I didn't have any values
            anomaly_e.to_netcdf(path = fileOUT_MERRA, mode ='w', engine='scipy')
            anomaly_e.close()
            
            #Create anomaly file
            anomaly_create = np.subtract(open_eto.ETo, anomaly_e.ETo_mean).rename('ETo_anom').to_dataset()
            anomaly_create.to_netcdf(fileOUT_MERRA_eto_anomaly, mode ='w', engine='scipy')
            anomaly_create.close()
            
    return(0)

#%% Windspeed

def run_MERRA_windspeed(wind):
    anomaly_r = xr.zeros_like(wind) #OUTPUT

    
    open_rzsm = wind #just rename
    try:
        check_mean_file=xr.open_dataset(fileOUT_MERRA_windspeed) #see if file is already created
        print(f'Already created MERRA windspeed anomaly file into {MERRA_dir}.')
    except FileNotFoundError:

        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for windspeed anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = fileOUT_MERRA_windspeed, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('windspeed_anom')
                anomaly_create.to_netcdf(fileOUT_MERRA_windspeed_anomaly, mode ='w', engine='scipy')
                anomaly_create.close()
                    

    return(0)
#%% Radiation
def run_MERRA_radiation(radiation):
    
    net_rad = np.add(radiation.LWGNET,radiation.SWNETSRF).to_dataset(name='srad')
    anomaly_r = xr.zeros_like(net_rad) #OUTPUT
    
    open_rzsm = net_rad#just rename
    try:
        check_mean_file=xr.open_dataset(fileOUT_MERRA_radiation) #see if file is already created
        print(f'Already created MERRA radiation anomaly file into {MERRA_dir}.')
    except FileNotFoundError:

        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for windspeed anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = fileOUT_MERRA_radiation, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('radiation_anom')
                anomaly_create.to_netcdf(fileOUT_MERRA_radiation_anomaly, mode ='w', engine='scipy')
                anomaly_create.close()
                
#%% Temperature max
def run_MERRA_tmax(temp):
    
    temp_ = temp.T2MMAX
    anomaly_r = xr.zeros_like(temp_).to_dataset() #OUTPUT
    
    open_rzsm = temp_.to_dataset()#just rename
    
    try:
        check_mean_file=xr.open_dataset(fileOUT_MERRA_tasmax) #see if file is already created
        print(f'Already created MERRA radiation anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for tasmax anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = fileOUT_MERRA_tasmax, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('tasmax_anom')
                anomaly_create.to_netcdf(fileOUT_MERRA_tasmax_anomaly, mode ='w', engine='scipy')
                anomaly_create.close()
#%% Temperature min
def run_MERRA_tmin(temp):
    f1 = fileOUT_MERRA_tasmin
    f2 = fileOUT_MERRA_tasmin_anomaly
    
    temp_ = temp.T2MMIN
    anomaly_r = xr.zeros_like(temp_).to_dataset() #OUTPUT
    
    open_rzsm = temp_.to_dataset()#just rename
    
    try:
        check_mean_file=xr.open_dataset(f1) #see if file is already created
        print(f'Already created MERRA radiation anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for tasmin anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = f1, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('tasmin_anom')
                anomaly_create.to_netcdf(f2, mode ='w', engine='scipy')
                anomaly_create.close()

#%% Temperature mean
def run_MERRA_tavg(temp):
    f1 = fileOUT_MERRA_tas
    f2 = fileOUT_MERRA_tas_anomaly
    
    temp_ = temp.T2MMEAN
    anomaly_r = xr.zeros_like(temp_).to_dataset() #OUTPUT
    
    open_rzsm = temp_.to_dataset()#just rename
    
    try:
        check_mean_file=xr.open_dataset(f1) #see if file is already created
        print(f'Already created MERRA radiation anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for tavg anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = f1, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('tas_anom')
                anomaly_create.to_netcdf(f2, mode ='w', engine='scipy')
                anomaly_create.close()

#%% Actual vapor pressure
def run_MERRA_avp(vapor_pressure):
    f1 = fileOUT_MERRA_actual_vapor_pressure
    f2 = fileOUT_MERRA_actual_vapor_pressure_anomaly
    
    avp_ = vapor_pressure.vapor_pressure
    anomaly_r = xr.zeros_like(avp_).to_dataset() #OUTPUT
    
    open_rzsm = avp_.to_dataset()#just rename
    
    try:
        check_mean_file=xr.open_dataset(f1) #see if file is already created
        print(f'Already created MERRA radiation anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
        for lat in range(open_rzsm.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_rzsm.Y.shape[0]))} for actual vapor pressure anomaly MERRA2.')
            for lon in range(open_rzsm.X.shape[0]):
                
                #Only work on areas within CONUS mask
                if CONUS_mask.CONUS_mask[0,lat,lon].values ==1:
                    
                    # all_values = open_f['ETo_gridmet'].isel(lat=lat,lon=lon)
                    #Choose the 7th day as a starting point because we need weekly data (start from 7, work backwards)
                    #Only need to look at 1 years worth of data because we are adding 
                    for day in range(7,377):
        
                        #Calculate anomaly based on a 42 day window on both sides of the 
                        #day (based on Julian day only)
                        #find the date of the current step
                        day_val=open_rzsm.time[day].to_numpy()
                  
                        yearly_average = {}
                        new_day_val = day_val
                        all_values = []
        
                        #Total of 22 years worth of data from SubX
                        for year_ in np.arange(2000,2023):
                            #Because we haven't taken a rolling mean yet, we need to here (only for anomaly - not day to day forecast comparisons)
                            
                            weekly_mean = np.nanmean(open_rzsm[name_1(open_rzsm)].isel(Y=lat,X=lon).sel(time=slice(day_val-np.timedelta64(7,'D'),day_val)))
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
                        day_val=open_rzsm.time[day].to_numpy()
                        for year_ in np.arange(2000,2023):
                            all_values.append(list(open_rzsm[name_1(open_rzsm)].sel(time=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(Y=lat,X=lon).to_numpy()))

                            if pd.to_datetime(day_val).year % 4 ==0:
                                day_val = day_val+np.timedelta64(366,'D')
                            else:
                                day_val = day_val+np.timedelta64(365,'D')
                            # open_f.ETo_gridmet.sel(day=slice(day_val-np.timedelta64(anomaly_range, 'D'),day_val+np.timedelta64(anomaly_range, 'D'))).isel(lat=lat,lon=lon).to_numpy()
                        
                        all_values.append(yearly_avg_list)
                        #flatten the list of lists, find mean
                        mean_ = np.nanmean([x for xs in all_values for x in xs])
                        
                        #Add the mean_ value to another dataset to use for the anomaly correlation coefficient
           
                        #Now add to the empty anomly_f file
                        for i in yearly_average.keys():
                            _date = pd.to_datetime(i)
                            #find the index in anomaly_date_list
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_r_list.index(_date)
                                anomaly_r[name_1(anomaly_r)][index_val, lat,lon] = mean_
                            except ValueError:
                                pass
                            
                
                anomaly_r[name_1(anomaly_r)][0:7,:,:] = anomaly_r[name_1(anomaly_r)][0+366:7+366].values #account for the first year because I didn't have any values
                anomaly_r.to_netcdf(path = f1, mode ='w', engine='scipy')
                anomaly_r.close()
                
                #Create anomaly file
                anomaly_create = np.subtract(open_rzsm[name_1(open_rzsm)], anomaly_r[name_1(anomaly_r)]).rename('avp_anom')
                anomaly_create.to_netcdf(f2, mode ='w', engine='scipy')
                anomaly_create.close()

#%%
if __name__ == '__main__':
    make_MERRA_ETo()
    run_MERRA_ETo_anomaly()
    run_MERRA_RZSM(open_rzsm)
    run_MERRA_windspeed(wind)
    run_MERRA_radiation(radiation)
    run_MERRA_tmax(temp)
    run_MERRA_tmin(temp)
    run_MERRA_tavg(temp)
    run_MERRA_avp(vapor_pressure)
    # run_GLEAM()
