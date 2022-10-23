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
import refet
from pyeto import fao


dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
MERRA_dir = f'{dir1}/Data/MERRA2'

evap='Penman'
eto_merra_save_name = f'{MERRA_dir}/ETo_{evap}_merged.nc'

fileOUT_MERRA = f'{MERRA_dir}/ETo_anomaly_{evap}_MERRA_mean.nc'
fileOUT_MERRA_eto_anomaly = f'{MERRA_dir}/ETo_anomaly_{evap}_MERRA.nc'

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

radiation = xr.open_dataset(f'{MERRA_dir}/radiation.nc4')
total_radiation = np.multiply((radiation.SWNETSRF + radiation.LWGNET).to_dataset(name='net_radiation'), 0.0864)  #convert from W/m2 to Mj/m2/d

radiation.close()
temp = np.subtract(xr.open_dataset(f'{MERRA_dir}/temperature.nc4'),273.15) #open, convert to Celsius


elevation = xr.open_dataset(f'{dir1}/Data/elevation/elev_regrid.nc',decode_times=False)

atm_pressure = (101.3 * ((293 - 0.0065 * elevation)/293)**5.26) # kPa --- LOOKS FINE AFTER INSPECTION


wind_humidity = xr.open_dataset(f'{MERRA_dir}/wind_humidity.nc4')

#QLML humidity MERRA2 is 
'''Calculate dewpoint from specific humidity'''
##' Convert specific humidity to relative humidity
##'
##' converting specific humidity into relative humidity
##' NCEP surface flux data does not have RH
##' from Bolton 1980 The computation of Equivalent Potential Temperature 
##' \url{http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html}
##' @title qair2rh
##' @param qair specific humidity, dimensionless (e.g. kg/kg) ratio of water mass / total air mass
##' @param temp degrees C
##' @param press pressure in mb
##' @return rh relative humidity, ratio of actual water mixing ratio to saturation mixing ratio
##' @export
##' @author David LeBauer
def qair2rh (qair, temp, press):
    press = np.multiply(press,10) #convert from kPa to millibars
    es =  6.112 * np.exp((17.67 * temp)/(temp + 243.5))
    e = qair * press / (0.378 * qair + 0.622)
    rh = e / es
    rh[rh > 1] <- 1
    rh[rh < 0] <- 0
    return(rh)

ea = xr.zeros_like(temp.T2MMAX)

#%%
def make_MERRA_ETo():
    try:
        save_output=xr.open_dataset(eto_merra_save_name)
    except FileNotFoundError:

        save_output = xr.zeros_like(total_radiation).rename(net_radiation='ETo')
        
        
        #need to only loop through days because ptC doesn't have any dates in the file
        for day in range(save_output.time.shape[0]):
            if day%100==0:
                print(f'Completed number {day} out of {max(range(save_output.time.shape[0]))}')

            for i_Y in range(save_output.Y.shape[0]):
                for i_X in range(save_output.X.shape[0]):
                    #Don't calculate extra data
                    if CONUS_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                        # print(i_X)
                        date_val = pd.to_datetime(save_output.time[day].values)
                        day_doy = date_val.timetuple().tm_yday #julian day
                        
                        rh = qair2rh(qair=wind_humidity.QLML[day, i_Y, i_X].values, \
                                     press=fao.atm_pressure(elevation.data[0,i_Y, i_X].values), \
                                         temp = temp.T2MMEAN[day, i_Y, i_X].values)
                                        

                        save_output.ETo[day, i_Y, i_X] = \
                                    fao.fao56_penman_monteith(net_rad = total_radiation.net_radiation[day, i_Y, i_X].values,\
                                                              t = temp.T2MMEAN[day, i_Y, i_X].values + 273.15, \
                                                              ws= wind_humidity.SPEED[day, i_Y, i_X].values, \
                                                              svp = fao.svp_from_t(temp.T2MMEAN[day, i_Y, i_X].values), \
                                                               avp = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(temp.T2MMIN[day, i_Y, i_X].values), \
                                                                                      svp_tmax=fao.svp_from_t(temp.T2MMAX[day, i_Y, i_X].values), \
                                                                                      rh_mean = rh),
                                                              delta_svp = fao.delta_svp(temp.T2MMEAN[day, i_Y, i_X].values), \
                                                              psy = fao.psy_const(fao.atm_pressure(elevation.data[0,i_Y, i_X].values)))
                                        
                                        
                        # print(fao.fao56_penman_monteith(net_rad = total_radiation.net_radiation[day, i_Y, i_X].values,\
                        #                                                               t = temp.T2MMEAN[day, i_Y, i_X].values + 273.15, \
                        #                                                               ws= wind_humidity.SPEED[day, i_Y, i_X].values, \
                        #                                                               svp = fao.svp_from_t(temp.T2MMEAN[day, i_Y, i_X].values), \
                        #                                                                avp = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(temp.T2MMIN[day, i_Y, i_X].values), \
                        #                                                                                       svp_tmax=fao.svp_from_t(temp.T2MMAX[day, i_Y, i_X].values), \
                        #                                                                                       rh_mean = rh),
                        #                                                               delta_svp = fao.delta_svp(temp.T2MMEAN[day, i_Y, i_X].values), \
                        #                                                               psy = fao.psy_const(fao.atm_pressure(elevation.data[0,i_Y, i_X].values))))
            
            
                        ea[day, i_Y, i_X] = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(temp.T2MMIN[day, i_Y, i_X].values), \
                                               svp_tmax=fao.svp_from_t(temp.T2MMAX[day, i_Y, i_X].values), \
                                                   
                                               rh_mean = rh)


        
        ea.to_dataset(name='vapor_pressure').to_netcdf(f'{MERRA_dir}/actual_vapor_pressure.nc4')

        save_output.to_netcdf(f'{eto_merra_save_name}')
        #Create a yearly sum file for inspection
        os.system(f'cdo monmean {eto_merra_save_name} {MERRA_dir}/eto_{evap}_monthly_mean.nc')
        os.system(f'cdo timmean -yearsum {eto_merra_save_name} {MERRA_dir}/eto_{evap}_mean_annual.nc')
        os.system(f'cdo -yearsum {eto_merra_save_name} {MERRA_dir}/eto_{evap}_year_average.nc')
        
        return(save_output)

#%%
#ETo anomalies
def run_MERRA_ETo_anomaly():
    open_eto = xr.open_dataset(f'{eto_merra_save_name}')
    #Create empty files rzsm
    anomaly_e = xr.zeros_like(open_eto)
    anomaly_e_mean = xr.zeros_like(anomaly_e)
    
    anomaly_date_e = pd.DataFrame(open_eto[list(open_eto.keys())[0]].time)
    anomaly_date_e.index = pd.to_datetime(anomaly_date_e.iloc[:,0])
    anomaly_date_e_list = list(anomaly_date_e.iloc[:,0])

    try:
        anomaly_e=xr.open_dataset(fileOUT_MERRA) #see if file is already created
        print(f'Already created MERRA anomaly file into {MERRA_dir}.')
    except FileNotFoundError:
    
        for lat in range(open_eto.Y.shape[0]):
            print(f'Working on lat {lat} out of {np.max(range(open_eto.Y.shape[0]))} for ETo anomaly MERRA2.')
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
                        day_val=open_eto.time[day].to_numpy()
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
                                index_val = anomaly_date_e_list.index(_date)
                                anomaly_e.ETo[index_val, lat,lon] = yearly_average[i]
                            except ValueError:
                                #If outside of array, don't process
                                pass
                            
                            try:
                                index_val = anomaly_date_e_list.index(_date)
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

#%%
if __name__ == '__main__':
    # run_MERRA_RZSM(anomaly_r)
    make_MERRA_ETo()
    run_MERRA_ETo_anomaly()
    # run_GLEAM()

