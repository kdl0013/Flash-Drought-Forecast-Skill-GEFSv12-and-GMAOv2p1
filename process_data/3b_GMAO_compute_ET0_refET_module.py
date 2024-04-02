#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate reference evapotranspiration (ETo) using ASCE Penman Monteith
    method https://www.apogeeinstruments.com/content/EWRI-ASCE-Reference-ET.pdf



@author: kdl
"""

#TODO: Working on finishing the ETo calculation,.... opening elevation file


import xarray as xr
import numpy as np
from numba import njit,prange
import os
import datetime as dt
import pandas as pd
from glob import glob
from numba import guvectorize, float64, int64


home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/lagged_ensemble_by_variable'

#TODO: Change code below to be a 'sed' from linux command
# dir1 = 'home_dir'
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/GMAO/lagged_ensemble_by_variable'
elevation_dir = f'{dir1}/Data/elevation/'
os.chdir(home_dir)


###Files
file_list = os.listdir()

soil_moisture = xr.open_dataset(glob('mrso*')[0]).stack(grid = ['X','Y'])
#Find nan locations for later use
nan_dataset = np.where(np.isnan(soil_moisture.READ_DESCRIPTION_ATTRS.values), np.nan, 0)

#Open file and change all values of 0.0 to nan
temperature_1 = (xr.open_dataset(glob('tas*')[0]).stack(grid = ['X','Y']))
# temperature_1 = temperature_1.where(temperature_1 != 0)

windV = (xr.open_dataset(glob('vas*')[0]).stack(grid = ['X','Y']))
# windV = windV.where(windV != 0)

windU = (xr.open_dataset(glob('uas*')[0]).stack(grid = ['X','Y']))
# windU = windU.where(windU != 0)

precip = (xr.open_dataset(glob('pr*')[0]).stack(grid = ['X','Y']))
# precip = precip.where(precip != 0)

dewpoint_temp_1 = (xr.open_dataset(glob('tdps*')[0]).stack(grid = ['X','Y']))
# dewpoint_temp_1= dewpoint_temp_1.where(dewpoint_temp_1 != 0)

shortwave_radiation = (xr.open_dataset(glob('dswrf*')[0]).stack(grid = ['X','Y']))
# shortwave_radiation = shortwave_radiation.where(shortwave_radiation != 0)

elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False).stack(grid = ['lon','lat'])
elevation = elevation.data[0,:]

dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values
#%%Temp is good, esat is good, delta looks good, 
#Step 1, Compute ET0 using the lagged average ensemble

#Constants
Cn = 1600 #mm/d
Cd = 0.38 #mm/s

def convert_temperature(temperature, dewpoint_temp) -> float:
    '''Convert temperature to Celsius (from Kelvin)'''
    temperature = np.subtract(temperature_1, 273.15)
    dewpoint_temp = np.subtract(dewpoint_temp_1, 273.15)
    
    return(temperature, dewpoint_temp)

temperature, dewpoint_temp = convert_temperature(temperature = temperature_1, \
                                                 dewpoint_temp = dewpoint_temp_1)

# #Check outputs
# print(temperature.READ_DESCRIPTION_ATTRS.values)
# print(temperature.READ_DESCRIPTION_ATTRS[0,200,0].values)
# print(np.nanmax(temperature.READ_DESCRIPTION_ATTRS.values))
# print(np.nanmin(temperature.READ_DESCRIPTION_ATTRS.values))

# print(dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values)
# print(dewpoint_temp.READ_DESCRIPTION_ATTRS[0,200,0].values)
# print(np.nanmax(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
# print(np.nanmin(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
##############################################################################
test_temp = 20
def esat_vapor_pressure(temperature, test_temp) -> float:
    '''Tetens Formula: Need temp in Celsius. Returns values in kPa
    Test temp should return 2339Pa or 2.34kPa'''
    esat = 0.6112 * np.exp((17.27*temperature)/(237.3 + temperature))
    esat_test = 0.6112 * np.exp((17.27*test_temp)/(237.3 + test_temp))
    delta = (2504 * esat) / (temperature + 237.3)**2
    return(esat, delta, esat_test)

esat, delta, esat_test = esat_vapor_pressure(temperature = temperature, test_temp=test_temp)

print()
print(f'Esat_test kPa is {esat_test} at {test_temp} C.')
print()
# print(esat.READ_DESCRIPTION_ATTRS.values)
print(f'Max esat value is {np.nanmax(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
print()
print(f'Min esat value is {np.nanmin(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
print()

# print(delta.READ_DESCRIPTION_ATTRS.values)
print(f'Max delta value is {np.nanmax(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
print()
print(f'Min delta value is {np.nanmin(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
print()
##############################################################################
test_dewpoint_temp = 16
def e_vapor_pressure(esat, dewpoint_temp, test_dewpoint_temp) -> float:
    '''Need temp in Celsius. Returns vapor pressure in Pa and humidity
    values in percent (e.g., 0.78)
    ea_test should return 1.81kPa'''
    ea = 0.6108 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
    ea_test = 0.6108 * np.exp((17.27*test_dewpoint_temp)/(237.3 + test_dewpoint_temp))
    rel_humidity = np.divide(ea,esat)
    return(ea, rel_humidity, ea_test)

ea, rel_humidity, ea_test = e_vapor_pressure(esat = esat, dewpoint_temp = dewpoint_temp,test_dewpoint_temp = 16)

print(f'Ea_test kPa is {ea_test} at {test_dewpoint_temp} C.')
print()
# print(ea.READ_DESCRIPTION_ATTRS.values)
print(f'Max ea value is {np.nanmax(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
print()
print(f'Min ea value is {np.nanmin(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
print()

# print(rel_humidity.READ_DESCRIPTION_ATTRS.values)
print(f'Max rel_humidity value is {np.nanmax(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
print()
print(f'Min rel_humidity value is {np.nanmin(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
print()
##############################################################################
def convert_dswrf_RN(shortwave_radiation) -> float:
    '''Returns shorwave radiation in MJ/m2 from W/m2'''
    shortwave_rad = shortwave_radiation / 1000000
    return(shortwave_rad)

shortwave_rad = convert_dswrf_RN(shortwave_radiation = shortwave_radiation)

# print(shortwave_rad.READ_DESCRIPTION_ATTRS.values)
print(f'Max shortwave_rad value is {np.nanmax(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
print()
print(f'Min shortwave_rad value is {np.nanmin(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
print()
##############################################################################
#Square both u and v and then take the square root
def compute_windSpeed(windU, windV) -> float:
    ''' Returns average daily wind at 10-m height m/s'''
    #Square u and v and take the square root to get windspeed
    windSpeed = xr.ufuncs.sqrt(np.square(windU) + np.square(windV))
    return (windSpeed)

windSpeed = compute_windSpeed(windU = windU, windV = windV)

# print(windSpeed.READ_DESCRIPTION_ATTRS.values)
print(f'Max windSpeed value is {np.nanmax(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
print()
print(f'Min windSpeed value is {np.nanmin(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
print()


##############################################################################
def mean_atmospheric_pressure(elevation) -> float:
    pressure = xr.zeros_like(temperature)

    for k in range(elevation.data.shape[0]):
        if k == 0 and i == 0 and j == 0:
            print(f'Pressure at elevation {elevation[k].values} meters is {101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)} kPa.')
        pressure.READ_DESCRIPTION_ATTRS[:,:,k] = 101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)
        # pressure = pressure.where(temperature_1 != 0)

    psychro = 0.000665 * pressure
    '''Returns atmospheric pressure in kPa. And returns psychrometric constant'''
    
    return(pressure, psychro)

pressure, psychro = mean_atmospheric_pressure(elevation = elevation)


# print(pressure.READ_DESCRIPTION_ATTRS.values)
print()
print(f'Max pressure value is {np.nanmax(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
print()
print(f'Min pressure value is {np.nanmin(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
print()

# print(psychro.READ_DESCRIPTION_ATTRS.values)
print(f'Max psychro value is {np.nanmax(psychro.READ_DESCRIPTION_ATTRS.values)}')
print()
print(f'Min psychro value is {np.nanmin(psychro.READ_DESCRIPTION_ATTRS.values)}')
print()

#%%
#TODO
'''NEED to convert all negative numbers to 0'''
def compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, \
                     esat, ea):
    
    # ETref = xr.zeros_like(windSpeed)
    '''Returns Reference Evapotranspiration in mm/d'''
    ETref = ((0.408 * delta * shortwave_radiation) + (psychro*(Cn/temperature+273) * windSpeed * \
        (esat - ea))) / (delta + (psychro * (1 + Cd * windSpeed)))
        
    return(ETref)
    
ETref_func = compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, esat, ea)

ETref_func.READ_DESCRIPTION_ATTRS.values

# print(ETref_func.READ_DESCRIPTION_ATTRS.values)
print()
print(f'Max ETref_func value is {np.nanmax(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
print()
print(f'Min ETref_func value is {np.nanmin(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
print()


ETref_unstack = ETref_func.unstack('grid').transpose('realization', 'time_series_dates', 'Y', 'X')

#Convert to an xarray object
var_OUT = xr.Dataset(
    data_vars = dict(
        ETo = (['realization','time_series_dates','Y','X'], ETref_unstack.READ_DESCRIPTION_ATTRS.values),
    ),
    coords = dict(
        X = ETref_unstack.X.values,
        Y = ETref_unstack.Y.values,
        time_series_dates = ETref_unstack.time_series_dates.values,
        realization = ETref_unstack.realization.values,
    ),
    attrs = dict(
        Description = 'Reference crop evapotranspiration (mm/day)'),
)

var_OUT.ETo.values

#Save as a netcdf for later processing
var_OUT.to_netcdf(path = f'{home_dir}/ETo_lagged_average.nc4', mode ='w')
#%%
a = ETref_func.READ_DESCRIPTION_ATTRS[0,200,0].values
print(a)


#%%


#%%

for j in range(test.READ_DESCRIPTION_ATTRS.shape[1]):
    for i in range(test.READ_DESCRIPTION_ATTRS.shape[2]):
        if test.READ_DESCRIPTION_ATTRS[0,j,i] < 0:
            print (test.READ_DESCRIPTION_ATTRS[0,i,i].values)
            break
        print (test.READ_DESCRIPTION_ATTRS[0,i,i].values) 
        
        
    



