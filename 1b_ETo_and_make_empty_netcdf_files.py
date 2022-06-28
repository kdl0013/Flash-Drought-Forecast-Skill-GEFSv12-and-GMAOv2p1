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

dir1 = 'main_dir'
num_processors = int('procs')
mod = 'model_name'
var = 'variable'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# num_processors = 10
# mod = 'GMAO'
# var = 'EDDI'

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
#All files have the same initialized days (part of the pre-processing that is 
#completed)
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc')):
        date_list.append(file[-13:-3])
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
        xr.open_dataset(f'ETo_{_date}.nc4')

        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
            
        print(f'Working on date {_date} to calculate SubX ETo and save into {home_dir}.')
        #Open up each file
        tasmax = xr.open_dataset(glob(f'tasmax*{_date}.nc')[0])
        tasmin = xr.open_dataset(glob(f'tasmin*{_date}.nc')[0])
        tdps = xr.open_dataset(glob(f'tdps*{_date}.nc')[0])
        dswrf = xr.open_dataset(glob(f'dswrf*{_date}.nc')[0])
        windU = xr.open_dataset(glob(f'uas*{_date}.nc')[0])
        windV = xr.open_dataset(glob(f'vas*{_date}.nc')[0])
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
                lead = output_f.L.values,
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
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)


#%%
'''Compute Reference ET for gridMET data.
 Because gridMET refET in its orginial state only contains 2 significant figures,
 i.e., 1.3 mm/d, we need to compute refET using gridMET data. But also, all gridmet
 data is 2 sig. fig. length, but this will create a float value.
 '''
def Refet_gridMET():
    
    try:
        xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
        # print('Already completed gridMET Reference ET.')
    except FileNotFoundError:
        
        print(f'Working on computing ET ref from gridMET variables and saving into {gridMET_dir}.')
        slice1 = "1999-01-10"
        slice2 = "2016-02-09"
        # gridMET_file = gridMET_file.sel(day=slice(f"{slice1}",f"{slice2}"))
        #Open up each file
        tasmax = xr.open_dataset(f'{gridMET_dir}/tmmx_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        tasmin = xr.open_dataset(f'{gridMET_dir}/tmmn_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        tavg = xr.open_dataset(f'{gridMET_dir}/tavg_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        dswrf = xr.open_dataset(f'{gridMET_dir}/srad_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #already in W/m2
        wind = xr.open_dataset(f'{gridMET_dir}/vs_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
        RHmin = xr.open_dataset(f'{gridMET_dir}/rmax_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
        RHmax = xr.open_dataset(f'{gridMET_dir}/rmin_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
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
        min_out, DN, DN = esat_vapor_pressure(tasmin, test_temp)
        max_out ,DN, DN = esat_vapor_pressure(tasmax, test_temp)
        
        ea = xr.zeros_like(tavg)
        ea = ((min_out.air_temperature * (RHmax.relative_humidity/100)) + (max_out.air_temperature * (RHmin.relative_humidity/100))) / 2
            

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
print('')
print('')
#Run function, saved
Refet_gridMET()
print('')
print('')
#%% Create empty EDDI.npy files for next script

#I brought this function out of loop, no need to repeat it more than once
T_FILE = xr.open_dataset(glob('ETo_1999-01-10.nc4')[0])

def make_empty_nc_files(init_date_list,T_FILE, var):
    print(f'Making empty .nc files for {var}.')
    # count=0
    
    # for var in ['EDDI','RZSM','ETo']:
    for _date in init_date_list:
 
        st_day =  pd.to_datetime(_date)
        st_day_2 = st_day + dt.timedelta(days=1)
         # day_doyey = st_day.timetuple().tm_yday #julian day
    
         #add more julian dates
        day_julian_a = [pd.to_datetime(st_day_2) + dt.timedelta(days=i) for i in range(45)]
        
        model_dirs = {}
        
        # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
        if var == 'EDDI':
            filename = 'EDDI'
            desc = 'Evaporative Demand Drought Index'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                
        elif var == 'RZSM':
            filename = f'{var}_anomaly'
            desc = f'RZSM anomaly SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
        elif var == 'ETo':
            filename = f'{var}_anomaly'
            desc = f'ETo anomaly SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
        
        # try:
        #     np.load(f'{home_dir}/julian_lead_{_date}.npy')
        # except FileNotFoundError:
        #     np.save(f'{home_dir}/julian_lead_{_date}.npy',day_julian_b) 
            
            
        #don't remake empty files
        try:
            xr.open_dataset(f'{home_dir}/{filename}_{_date}.nc4')
            
        except FileNotFoundError:
            #Make empty files to insert 
            empty = (np.zeros_like(T_FILE.to_array()).squeeze())         
            template = xr.open_dataset(f'ETo_{_date}.nc4')
    
            #initialize empty file
            var_OUT = xr.zeros_like(template)

            if var == 'RZSM':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        RZSM_anom = (['S','model','lead','Y','X'], empty[:,:,:,:,:]),
                    ),
                    coords = dict(
                        X = var_OUT.X.values,
                        Y = var_OUT.Y.values,
                        lead = var_OUT.lead.values,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                
                var_final.RZSM_anom[1,:,:,:,:] = np.nan
                #Drop S dimension to save storage space
                var_final = var_final.dropna(dim='S',how='all')
                var_final.RZSM_anom[0,:,:,:,:] = np.nan
                
                
            elif var == 'EDDI':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        EDDI = (['S','model','lead','Y','X'], empty[:,:,:,:,:]),
                    ),
                    coords = dict(
                        X = var_OUT.X.values,
                        Y = var_OUT.Y.values,
                        lead = var_OUT.lead.values,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                
                
                var_final.EDDI[1,:,:,:,:] = np.nan
                #Drop S dimension to save storage space
                var_final = var_final.dropna(dim='S',how='all')
                
            if var == 'ETo':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        ETo_anom = (['S','model','lead','Y','X'], empty[:,:,:,:,:]),
                    ),
                    coords = dict(
                        X = var_OUT.X.values,
                        Y = var_OUT.Y.values,
                        lead = var_OUT.lead.values,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                
                var_final.ETo_anom[1,:,:,:,:] = np.nan
                #Drop S dimension to save storage space
                var_final = var_final.dropna(dim='S',how='all')
                var_final.ETo_anom[0,:,:,:,:] = np.nan

            var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc')
            #compress so that when I re-write the file, it is quicker
            # '''But this doesn't work after the file is re-read and re-saved'''
            # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
            # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
            var_final.close()

   
#Create EDDI files           
make_empty_nc_files(init_date_list = init_date_list,T_FILE = T_FILE, var=var)
#%%'''Make empty files to keep track of anomaly calculations'''

#Make a completed list file for EDDI, add new names to next code block
new_eddi = f'EDDI_completed_nc_{mod}.txt'
new_eto_anom = f'ETo_completed_anomaly_nc_{mod}.txt'
new_rzsm = f'RZSM_completed_anomaly_nc_{mod}.txt'

#Create a new file for each index to keep track of what has been completed 
#since this is a pretty slow process
for name in [new_eddi,new_eto_anom,new_rzsm]:
    try:
        completed_dates = np.loadtxt(f'{script_dir}/{name}',dtype='str')         
    except OSError:
        os.system(f'touch {script_dir}/{name}')
        name_index = f'{name}'.split('_')[0]
        os.system(f'echo "Completed {name_index}" > {script_dir}/{name}')
