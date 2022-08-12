#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate reference evapotranspiration (ETo) using ASCE Penman Monteith
    method https://pypi.org/project/refet/ python package : refet
    
    gridMET calculate reference evapotranspiration with the same python package
    
    Calculate EDDI from SubX data (EDDI is already downloaded for historical data)
    


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
from datetime import timedelta
from metpy import calc

dir1 = 'main_dir'
num_processors = int('procs')
mod = 'model_name'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# num_processors = 10
# mod = 'GMAO'
# var = 'ETo'

home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'
gridMET_dir = f'{dir1}/Data/gridMET'
smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'

###Files
file_list = os.listdir()

#variables = ['dswrf','tasmax', 'tasmin', 'uas', 'vas']

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of _datethe pre-processing that is 
#completed)
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'vas')    
                  #max     #min   #dew   #rad   #wind
#open files [need tasmax, tasmin, tdps, dswrf, windSpeed]

# for _date in init_date_list:
#_date = init_date_list[0]
#%%
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX(_date):
    
    try:
        xr.open_dataset(f'ETo_{_date}.nc')

        print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
            
        print(f'Working on date {_date} to calculate SubX ETo and save into {home_dir}.')
        #Open up each file
        tasmax = xr.open_dataset(glob(f'tasmax*{_date}.nc4')[0])
        tasmin = xr.open_dataset(glob(f'tasmin*{_date}.nc4')[0])
        tdps = xr.open_dataset(glob(f'tdps*{_date}.nc4')[0])
        dswrf = xr.open_dataset(glob(f'dswrf*{_date}.nc4')[0])
        windU = xr.open_dataset(glob(f'uas*{_date}.nc4')[0])
        windV = xr.open_dataset(glob(f'vas*{_date}.nc4')[0])
        elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False)
    
        #CONUS mask
        conus_mask = xr.open_dataset(f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4')
        
    
        def convert_temperature(tasmax, tasmin, dewpoint_temp) -> float:
            '''Convert temperature to Celsius (from Kelvin), calculate vapor pressure (ea)
            from dewpoint temp'''
            tasmax = np.subtract(tasmax, 273.15)
            tasmin = np.subtract(tasmin, 273.15)
            dewpoint_temp = np.subtract(dewpoint_temp, 273.15)
            ea = 0.6112 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
            
            return(tasmax, tasmin, dewpoint_temp, ea)
        
        tasmax, tasmin, dewpoint_temp, ea = convert_temperature(tasmax = tasmax, tasmin = tasmin,
                                                         dewpoint_temp = tdps)     
    
        def convert_dswrf_RN(shortwave_radiation) -> float:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            shortwave_rad = shortwave_radiation / 1000000
            return(shortwave_rad)
    
        dswrf = convert_dswrf_RN(shortwave_radiation = dswrf)
    
        def compute_windSpeed(windU, windV) -> float:
            ''' Returns average daily wind at 10-m height m/s'''
            #Square u and v and take the square root to get windspeed
            windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
            return (windSpeed)
    
        windSpeed = compute_windSpeed(windU = windU, windV = windV)
        
        
        #convert to lead to julian day for later processing anomalies
        def date_file_info(SubX_file):

            a_date_in= SubX_file.L.values
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

            return(a_date_out)
        
        julian_list = date_file_info(tasmax)
        
        output_f = xr.zeros_like(tasmax)
        
        for i_mod in range(tasmax.tasmax.shape[1]):
            for i_lead in range(tasmax.tasmax.shape[2]):
                #changed these 2 lines of code from under i_X
  
                for i_Y in range(tasmax.tasmax.shape[3]):
                    for i_X in range(tasmax.tasmax.shape[4]):
                        
                        #Don't calculate extra data
                        if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:

                            # print(i_X)
                            date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                            day_doy = date_val.timetuple().tm_yday #julian day
                        
                            output_f.tasmax[0,i_mod, i_lead, i_Y, i_X] = \
                                (refet.Daily(tmin = tasmin.tasmin[0,i_mod, i_lead, i_Y, i_X].values, \
                                            tmax = tasmax.tasmax[0,i_mod, i_lead, i_Y, i_X].values, \
                                            ea = ea.tdps[0,i_mod, i_lead, i_Y, i_X].values, \
                                            rs = dswrf.dswrf[0,i_mod, i_lead, i_Y, i_X].values, \
                                            uz = windSpeed[0,i_mod, i_lead, i_Y, i_X].values, \
                                            zw = 10, elev = elevation.data[0,i_Y, i_X].values,
                                            lat = tasmax.Y[i_Y].values,
                                            doy = day_doy, method = 'asce',
                                            input_units={'tmin': 'C', 'tmax': 'C', \
                                                         'rs': 'Langleys', 'uz': 'm/s', \
                                                             'lat': 'deg'}).etr())[0]   
                            # print((refet.Daily(tmin = tasmin.tasmin[0,i_mod, i_lead, i_Y, i_X].values, \
                            #             tmax = tasmax.tasmax[0,i_mod, i_lead, i_Y, i_X].values, \
                            #             ea = ea.tdps[0,i_mod, i_lead, i_Y, i_X].values, \
                            #             rs = dswrf.dswrf[0,i_mod, i_lead, i_Y, i_X].values, \
                            #             uz = windSpeed[0,i_mod, i_lead, i_Y, i_X].values, \
                            #             zw = 10, elev = elevation.data[0,i_Y, i_X].values,
                            #             lat = tasmax.Y[i_Y].values,
                            #             doy = day_doy, method = 'asce',
                            #             input_units={'tmin': 'C', 'tmax': 'C', \
                            #                          'rs': 'Langleys', 'uz': 'm/s', \
                            #                              'lat': 'deg'}).etr())[0]   )
        
        #Sometime infinity values mess up the saving of the file
        output_f = output_f.where(output_f.apply(np.isfinite)).fillna(0.0)
        output_f.tasmax[1,:,:,:,:] = np.nan
        #Drop S dimension to save storage space
        output_f = output_f.dropna(dim='S',how='all')
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo = (['S', 'model','lead','Y','X'], output_f.tasmax.values),
            ),
            coords = dict(
                X = output_f.X.values,
                Y = output_f.Y.values,
                lead = julian_list,
                model = output_f.M.values,
                S = output_f.S.values
            ),
            attrs = dict(
                Description = 'Reference crop evapotranspiration (mm/day)'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{home_dir}/ETo_{_date}.nc', mode ='w')
        print(f'Saved {_date} into {home_dir}.')


#%%
'''Compute Reference ET for gridMET data.

 '''
def Refet_gridMET():
    
    try:
        xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
    except FileNotFoundError:
        
        print(f'Working on computing ET ref from gridMET variables and saving into {gridMET_dir}.')
        # gridMET_file = gridMET_file.sel(day=slice(f"{slice1}",f"{slice2}"))
        #Open up each file
        tasmax = xr.open_dataset(f'{gridMET_dir}/tmmx_remap_final.nc')#Kelvin
        tasmin = xr.open_dataset(f'{gridMET_dir}/tmmn_remap_final.nc')#Kelvin
        tavg = xr.open_dataset(f'{gridMET_dir}/tavg_remap_final.nc')#Kelvin
        dswrf = xr.open_dataset(f'{gridMET_dir}/srad_remap_final.nc')#already in W/m2
        wind = xr.open_dataset(f'{gridMET_dir}/vs_remap_final.nc')#m/s2, 10 ft.
        RHmin = xr.open_dataset(f'{gridMET_dir}/rmax_remap_final.nc')#m/s2, 10 ft.
        RHmax = xr.open_dataset(f'{gridMET_dir}/rmin_remap_final.nc')#m/s2, 10 ft.
        # gridMET_file = gridMET_file.sel(day=slice("1999-01-11","2016-02-09"))
            
        elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False)
        
        #CONUS mask
        conus_mask = xr.open_dataset(f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4')
        
        def convert_temperature(tasmax, tasmin) -> float:
            '''Convert temperature to Celsius (from Kelvin), calculate vapor pressure (ea)
            from dewpoint temp'''
            tasmax = np.subtract(tasmax, 273.15)
            tasmin = np.subtract(tasmin, 273.15)
            
            return(tasmax, tasmin)
        #Call function for temperature files
        tasmax, tasmin = convert_temperature(tasmax = tasmax, tasmin = tasmin)
                
        test_temp = 20
        def esat_vapor_pressure(tavg, test_temp) -> float:
            '''Tetens Formula: Need temp in Celsius. Returns values in kPa
            Test temp should return 2339Pa or 2.34kPa'''
            esat = 0.6112 * np.exp((17.27*tavg)/(237.3 + tavg))
            esat_test = 0.6112 * np.exp((17.27*test_temp)/(237.3 + test_temp))
            delta = (2504 * esat) / (tavg + 237.3)**2
            return(esat, delta, esat_test)
        
        #Calculate saturated vapor pressure:
        esat, delta, esat_test = esat_vapor_pressure(tavg, test_temp)
        #Calculate vapor pressure:
        esat_min, DN, DN = esat_vapor_pressure(tasmin, test_temp)
        esat_max ,DN, DN = esat_vapor_pressure(tasmax, test_temp)
        
        ea = xr.zeros_like(tavg)
        #Take the average of the actual vapor pressure
        #Source https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
        ea = ((esat_min.air_temperature * (RHmax.relative_humidity/100)) + \
              (esat_max.air_temperature * (RHmin.relative_humidity/100))) / 2
            

        def convert_dswrf_RN(shortwave_radiation) -> float:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            shortwave_rad = shortwave_radiation / 1000000
            return(shortwave_rad)
    
        dswrf = convert_dswrf_RN(shortwave_radiation = dswrf)
    
        output_f = xr.zeros_like(tasmax)

        for day_ in range(tasmax.day.shape[0]):
            for lat in range(tasmax.Y.shape[0]):
                for lon in range(tasmax.X.shape[0]):
                    day_doy = pd.to_datetime(tasmax.day[day_].values).timetuple().tm_yday #julian day
                    
                    output_f.air_temperature[day_, lat, lon] = \
                        (refet.Daily(tmin = tasmin.air_temperature[day_, lat, lon].values, \
                                    tmax = tasmax.air_temperature[day_, lat, lon].values, \
                                    ea = ea[day_, lat, lon].values, \
                                    rs = dswrf.surface_downwelling_shortwave_flux_in_air[day_, lat, lon].values, \
                                    uz = wind.wind_speed[day_, lat, lon].values, \
                                    zw = 10, elev = elevation.data[0,lat, lon].values,
                                    lat = tasmax.Y[lat].values,
                                    doy = day_doy, method = 'asce',
                                    input_units={'tmin': 'C', 'tmax': 'C', \
                                                 'rs': 'Langleys', 'uz': 'm/s', \
                                                     'lat': 'deg'}).etr())[0]                            
    
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_gridmet = (['day', 'lat','lon'], output_f.air_temperature.values),
            ),
            coords = dict(
                day = output_f.day.values,
                lat = output_f.Y.values,
                lon = output_f.X.values
            ),
            attrs = dict(
                Description = 'Reference crop evapotranspiration (mm/day)'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{gridMET_dir}/ETo_gridMET_merged.nc', mode ='w')
        print(f'Saved file into {gridMET_dir}.')
#%%    



# os.getcwd()
# #convert to lead to julian day for later processing anomalies
# def date_file_info(SubX_file):

#     a_date_in= SubX_file.lead.values
#     #get the start date
#     a_start_date = pd.to_datetime(SubX_file.S.values[0])
#     a_date_out=[]
#     for a_i in range(len(a_date_in)):
#         a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

#     return(a_date_out)
  

# for file in sorted(glob('ETo*.nc4')):
#     open_f = xr.open_dataset(file)
#     open_f.close()
#     julian_list = date_file_info(open_f)
#     open_f = open_f.assign_coords(lead=julian_list)
#     open_f.to_netcdf(f'{file}5')

#Save relative humidity into it's own dataset for each variable

def relative_humidity_subx(_date):

    try:
        xr.open_dataset(f'RelativeHumidity_{_date}.nc4')

        print(f'{_date} already completed for RH. Saved in {home_dir}.')
    except FileNotFoundError:
            
        print(f'Working on date {_date} to calculate SubX ETo and save into {home_dir}.')
        #Open up each file
        tasmax = xr.open_dataset(glob(f'tasmax*{_date}.nc4')[0])
        tasmin = xr.open_dataset(glob(f'tasmin*{_date}.nc4')[0])
        tdps = xr.open_dataset(glob(f'tdps*{_date}.nc4')[0])
    
# tasmax.tasmax[:,:,0:7,:,:].mean(dim='M').shape
        def convert_temperature(tasmax, tasmin, dewpoint_temp) -> float:
            '''Convert temperature to Celsius (from Kelvin), calculate vapor pressure (ea)
            from dewpoint temp'''
            tasmax = np.subtract(tasmax, 273.15)
            tasmin = np.subtract(tasmin, 273.15)
            dewpoint_temp = np.subtract(dewpoint_temp, 273.15)
            ea = 0.6112 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
            
            return(tasmax, tasmin, dewpoint_temp, ea)
        
        tasmax, tasmin, dewpoint_temp, ea = convert_temperature(tasmax = tasmax, tasmin = tasmin,
                                                         dewpoint_temp = tdps)     
        
        tavg = (tasmax.tasmax + tasmin.tasmin)/2
        
        #formula https://www.omnicalculator.com/physics/relative-humidity
        rh = 100 * (ea/(0.6112*(np.exp((17.27*tavg)/(237.3 + tavg)))))
            
        # rh.tdps.values
        #save to a new file
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                RH = (['S', 'model','lead','Y','X'], rh.tdps.values),
            ),
            coords = dict(
                X = tasmax.X.values,
                Y = tasmax.Y.values,
                lead = tasmax.L.values,
                model = tasmax.M.values,
                S = tasmax.S.values
            ),
            attrs = dict(
                Description = 'Relative humidity computed from temperature and dewpoint.'),
        )                    
    
        var_OUT.to_netcdf(path = f'{home_dir}/RelativeHumidity_{_date}.nc4', mode ='w')
        print(f'Saved {_date} into {home_dir}.')
    
#%%
'''Now we save relative humidity for gridMET'''
def relative_humidity_gridMET():
    
    try:
        xr.open_dataset(f'{gridMET_dir}/RelativeHumidity_gridMET_merged.nc')
    except FileNotFoundError:
        
        print(f'Saving relative humidity gridMET variables and saving into {gridMET_dir}.')

        # gridMET_file = gridMET_file.sel(day=slice(f"{slice1}",f"{slice2}"))
        #Open up each file
        RHmin = xr.open_dataset(f'{gridMET_dir}/rmax_remap_final.nc')#m/s2, 10 ft.
        RHmax = xr.open_dataset(f'{gridMET_dir}/rmin_remap_final.nc')#m/s2, 10 ft.
        # gridMET_file = gridMET_file.sel(day=slice("1999-01-11","2016-02-09"))
            
        #take the mean of RHmin and RHmax
        rh_final = (RHmin.relative_humidity + RHmax.relative_humidity)/2

        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                RH = (['day', 'lat','lon'], rh_final.values),
            ),
            coords = dict(
                day = RHmin.day.values,
                lat = RHmin.Y.values,
                lon = RHmin.X.values
            ),
            attrs = dict(
                Description = 'Relative Humidity from RHmin and RHmax'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{gridMET_dir}/RelativeHumidity_gridMET_merged.nc', mode ='w')
        print(f'Saved file into {gridMET_dir}.')
        
        
#%%    
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)
    print('')
    print('')
    #Run function, saved
    Refet_gridMET()
    print('')
    print('')
    p.map(relative_humidity_subx, init_date_list)
    print('')
    print('')
    relative_humidity_gridMET()
