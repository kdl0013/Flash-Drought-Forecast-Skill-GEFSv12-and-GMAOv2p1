#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    
    https://wetlandscapes.github.io/blog/blog/penman-monteith-and-priestley-taylor/
    
    https://soilwater.github.io/pynotes-agriscience/notebooks/evapotranspiration.html
    
    ## Priestley-Taylor (1972)
    change altitude to elevation
    change alpha to be priestely taylor coefficient

    def priestley_taylor(T_min,T_max,solar_rad,altitude):
        T_avg = (T_min + T_max)/2
        atm_pressure = 101.3 * ((293 - 0.0065 * altitude)/293)**5.26 # kPa
        gamma = 0.000665 * atm_pressure # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
        delta = 4098 * (0.6108 * np.exp(17.27 * T_avg / (T_avg  + 237.3))) / (T_avg  + 237.3)**2 # Slope saturated vapor pressure
        soil_heat_flux = 0;
        alpha = 0.5
        PET = alpha*delta/(delta+gamma)*(solar_rad-soil_heat_flux)
        return PET
        
    

    


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

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = int('7')
mod = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# num_processors = 10
# mod = 'GMAO'
# var = 'ETo'

home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'
ptC_dir = f'{dir1}/Data/Priestley_Taylor_Makkink_evap_coeff_maps'
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
        
init_date_list = return_date_list(var = 'tas')    
                  #max     #min   #dew   #rad   #wind
#open files [need tasmax, tasmin, tdps, dswrf, windSpeed]

# for _date in init_date_list:
# _date = init_date_list[0]
# 

'''IMPORTANT NOTE: RSMAS has no dlwrf data for year 2018 or 2019. So if no data, 
just use the shortwave flux alone as the net radiation'''

#%%
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX(_date):
# for _date in init_date_list:
    try:
        xr.open_dataset(f'ETo_{_date}.nc4')

        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
            
        print(f'Working on date {mod} {_date} to calculate Priestley-Taylor ETo.')
        #Open up each file
        try:
            tavg = xr.open_dataset(glob(f'tas*{_date}.nc4')[0])
            tavg = (np.subtract(tavg,273.15))  #convert to Celsius
        except ValueError:
            exit()
        
        #Missing data
        try:
            dswrf = xr.open_dataset(glob(f'dswrf*{_date}.nc4')[0])
            dswrf = np.divide(dswrf,1000000) #convert to MJ/m2
        except IndexError:
            pass

        try:
            dlwrf = xr.open_dataset(glob(f'dlwrf*{_date}.nc4')[0])
            dlwrf = np.divide(dlwrf,1000000) #convert to MJ/m2
        except IndexError:
            pass

        #Check if object is already in memory
        if ('dlwrf' in list(locals().keys())) and ('dswrf' in list(locals().keys())):
            srad = np.add(dswrf.dswrf,dlwrf.dlwrf).to_dataset()
        elif 'dswrf' in list(locals().keys()):
            srad = dswrf
        elif 'dlwrf' in list(locals().keys()):
            srad = dlwrf
        
        elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False)
        ptC = xr.open_dataset(f'{ptC_dir}/pt_coeff_FINAL.nc')

        #convert to lead to julian day for later processing anomalies
        def date_file_info(SubX_file):

            a_date_in= SubX_file.L.values
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

            return(a_date_out)
        
        julian_list = date_file_info(tavg)
        
        output_f = xr.zeros_like(tavg)
        
        def priestley_taylor(tavg1,srad1,elevation1,ptC1):
            atm_pressure = 101.3 * ((293 - 0.0065 * elevation1)/293)**5.26 # kPa
            gamma = 0.000665 * atm_pressure # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
            delta = 4098 * (0.6108 * np.exp(17.27 * tavg1 / (tavg1  + 237.3))) / (tavg1  + 237.3)**2 # Slope saturated vapor pressure
            soil_heat_flux = np.zeros_like(delta)
            PET = ptC1*delta/(delta+gamma)*(srad1-soil_heat_flux)
            return PET

        for i_mod in range(tavg.tas.shape[1]):
            for i_lead in range(tavg.tas.shape[2]):

                tavg1=tavg[list(tavg.keys())[0]].isel(L=i_lead,M=i_mod).values[0]
                srad1=srad[list(srad.keys())[0]].isel(L=i_lead,M=i_mod).values[0]
                elevation1=elevation[list(elevation.keys())[0]].isel(time=0).values
                ptC1=ptC[list(ptC.keys())[0]].values
                
                output_f.tas[0,i_mod, i_lead, :, :] = priestley_taylor(tavg1,srad1,elevation1,ptC1)*1000

        #Sometime infinity values mess up the saving of the file
        output_f = output_f.where(output_f.apply(np.isfinite)).fillna(0.0)
        
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo = (['S', 'model','lead','Y','X'], output_f.tas.values),
            ),
            coords = dict(
                X = output_f.X.values,
                Y = output_f.Y.values,
                lead = julian_list,
                model = output_f.M.values,
                S = output_f.S.values
            ),
            attrs = dict(
                Description = 'Reference crop evapotranspiration (mm/day). Priestley-Taylor formula'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{home_dir}/ETo_{_date}.nc4', mode ='w')
        print(f'Saved {_date} into {home_dir}.')


#%%
'''Compute Reference ET for gridMET data.

 '''
def Refet_MERRA():
    
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
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)
