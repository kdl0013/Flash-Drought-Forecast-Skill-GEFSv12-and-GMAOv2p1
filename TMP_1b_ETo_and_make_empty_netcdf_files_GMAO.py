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

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = int('10')
mod = 'GMAO'
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

        print(f'{_date} already completed for ETo. Saved in {home_dir}.')
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
            windSpeed = xr.ufuncs.sqrt((np.square(windU.uas) + np.square(windV.vas)))
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
        var_OUT.to_netcdf(path = f'{home_dir}/ETo_{_date}.nc4', mode ='w')
        print(f'Saved file into {home_dir}.')
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
        print('Already completed gridMET Reference ET.')
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

def make_empty_nc_files(init_date_list,T_FILE):
    print('Making empty .npy and .nc4 files for EDDI, RZSM, and ETo.')
    # count=0
    
    for var in ['EDDI','RZSM','ETo']:
        for _date in init_date_list:
     
            st_day =  pd.to_datetime(_date)
            st_day_2 = st_day + dt.timedelta(days=1)
             # day_doyey = st_day.timetuple().tm_yday #julian day
        
             #add more julian dates
            day_julian_a = [pd.to_datetime(st_day_2) + dt.timedelta(days=i) for i in range(45)]
            day_julian_b = [i.timetuple().tm_yday for i in day_julian_a]                            
            
            model_dirs = {}
            
            # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
            if var == 'EDDI':
                filename = 'EDDI'
                desc = 'Evaporative Demand Drought Index'
                for model in [0,1,2,3]:
                    model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                    
            elif var == 'RZSM':
                filename = f'{var}_anomaly'
                desc = 'RZSM anomaly. Calculated from 7-day average mean by julian day and lead week.'
                for model in [0,1,2,3]:
                    model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
                
            elif var == 'ETo':
                filename = f'{var}_anomaly'
                desc = 'ETo anomaly. Calculated from 7-day average mean by julian day and lead week.'
                for model in [0,1,2,3]:
                    model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
            try:
                np.load(f'{home_dir}/julian_lead_{_date}.npy')
            except FileNotFoundError:
                np.save(f'{home_dir}/julian_lead_{_date}.npy',day_julian_b) 
                
                
            #don't remake empty files
            try:
                xr.open_dataset(f'{home_dir}/{filename}_{_date}.nc4')
                
            except FileNotFoundError:
                #Make empty files to insert 
                empty = (np.zeros_like(T_FILE.to_array()).squeeze())
                # #Only do 1 model at a time
                # empty = empty[:,0,:,:,:]
                # empty[:,:,:,:] = np.nan
                # np.save(f'{home_dir}/EDDI_{_date}.npy',empty)
                # np.save(f'{home_dir}/ETo_anomaly_{_date}.npy',empty)
                # np.save(f'{home_dir}/RZSM_anomaly_{_date}.npy',empty)
                
                #save lead values as julian day for later processing

                         
                template = xr.open_dataset(f'ETo_{_date}.nc4')
        
                #initialize empty file
                var_OUT = xr.zeros_like(template)
                
                # #Add data from each model into the empty file
                # for model in model_dirs.items():
                #     var_OUT.ETo[:,int(model[0][-1]),:,:,:] = np.load(f'{filename}_{_date}.npy',allow_pickle=True)
                #     #print(mod)
                  
                
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

                var_final.to_netcdf(path = f'{home_dir}/{filename}_{_date}.nc4', mode ='w',engine='scipy')
                #compress so that when I re-write the file, it is quicker
                '''But this doesn't work after the file is re-read and re-saved'''
                os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
                os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
                
                var_final.close()

   
#Create EDDI files           
make_empty_nc_files(init_date_list = init_date_list,T_FILE = T_FILE)
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

#%%
'''Previous attempts to calculate reference ET'''


# #%%
# '''process SubX files and create EDDI values'''
# def multiProcess_EDDI_SubX_TEST(_date):
# # for _date in init_date_list[start_:end_]:
#     try:
#         #don't re-work code if the netcdf file is already created
#         xr.open_dataset(f'{home_dir}/EDDI_{_date}.nc4')
#     except FileNotFoundError:
#         print(f'Calculating EDDI on SubX for {_date} with all models and saving as .npy in {home_dir}.')
# # for _date in init_date_list[0:80]:   
#         os.chdir(f'{home_dir}')
#         week_lead = 6

#         #Used for eliminating iterating over grid cells that don't matter
#         smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
      
    
#         #get julian day, timestamps, and the datetime
#         def date_file_info(_date, variable):
#             open_f = xr.open_dataset(f'{variable}_{_date}.nc4')
#             a_date_in= open_f[f'{variable}'].lead.values
#             #get the start date
#             a_start_date = pd.to_datetime(open_f.S.values[0])
#             #Add the dates based on index value in date_in
#             a_date_out = []
    
#             for a_i in range(len(a_date_in)):
#                 a_date_out.append(a_start_date + dt.timedelta(days = a_i))
                
#             start_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
#             #Julian day into a list                            
#             a_julian_out = [start_julian + i for i in range(len(a_date_out))]
            
#             #month of file
#             INdate_for_month = dt.datetime(int(_date[0:4]),int(_date[5:7]),int(_date[8:10]))
            
#             return(open_f,a_date_out,a_julian_out,INdate_for_month)
        
#         subx,file_timestamp_list,file_julian_list,file_datetime = date_file_info(_date=_date,variable='ETo')
         
#         #Now convert to julian date and append coordinates
#         subx2 = subx.assign_coords(lead = file_julian_list)

#         for i_Y in range(subx2.ETo.shape[3]):
#             for i_X in range(subx2.ETo.shape[4]):
#                 if _date == '1999-01-10':
#                     print(f'Working on lat {i_Y} and lon {i_X}')

                              
#                 '''Next if statement creates a mask to restrict within CONUS, because
#                 EDDI covers North America, this will create the mask using SMERGE plus 2 additional grid cell.
#                 Don't loop over the values that shouldn't have values.
#                 SMERGE indexes 38,6 and 38,7 are missing, but we can compute EDDI values '''
#                 if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

#                     def dict1_subx2():
#                         #append 7-day summation from all files to a new dictionary, 
#                         #trying to avoid having to loop through models.
#                         summation_ETo_mod0 = {}
#                         summation_ETo_mod1 = {}
#                         summation_ETo_mod2 = {}
#                         summation_ETo_mod3 = {}
#                         models = [0,1,2,3]
                        
#                         for julian_d in file_julian_list:
#                             #You must julian_d + week_lead because with EDDI you need 7-day looking backwards into the past. Must have 7 values.
#                             #Choose just model = 0 because we just need to know if there are 7-days total in any model
#                             if (len(subx2.ETo.sel(lead=slice(julian_d,week_lead + julian_d)).isel(S=0, model=0, X=i_X, Y=i_Y).values)) == 7:
#                                 summation_ETo_mod0[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod1[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod2[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod3[f'{week_lead+julian_d}']=[]

#                                 summation_ETo_mod0[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[0], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod1[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[1], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod2[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[2], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod3[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[3], X=i_X, Y=i_Y).values)})   

#                         '''7-day sum of Eto by index:
#                         Next we will append to each julian day value in the key in the dictionary with
#                         the same julian day from all files'''
                        
#                         #Keep only files within a certain date
#                         def limit_file_list(_date):
#                             month_start = pd.to_datetime(_date).month
                            
#                             dates_to_keep = []
#                             for file in sorted(glob('ETo*.nc4')):
#                                 test_1 = month_start - pd.to_datetime(file[-14:-4]).month
#                                 test_2 = pd.to_datetime(file[-14:-4]).month - month_start
                                
#                                 if test_1 <=-9 or test_1>= 9:
#                                     dates_to_keep.append(file)
#                                 elif abs(test_1) <=3 or abs(test_2) <=3:
#                                     dates_to_keep.append(file)
                                    
#                             return (dates_to_keep)
                        
#                         dates_to_keep = limit_file_list(_date)
                        
#                         for file in dates_to_keep:
#                             #Remove unnecessary files that won't contain any of the dates, if difference in months
#                             #is >=2, then skip that file because of 45 day lead time.
                            
#                             # OUTdate_for_month = dt.datetime(int(file[4:8]),int(file[9:11]),int(file[12:14]))
                                                        
#                             # num_months = (OUTdate_for_month.year - file_datetime.year) + (OUTdate_for_month.month - file_datetime.month) 
#                             # (OUTdate_for_month.month - file_datetime.month)
                            
#                             #Dont' re-open the same file
#                             if file[-14:-4] != _date:
                                
#                                 open_f = xr.open_dataset(file)
                                
#                                 '''Convert lead dates to a vector, then add it back into a netcdf
#                                 because I cannot convert np.datetime to a pd.datetime, I will need
#                                 to iterate over each of the dates in a list from Et_ref'''
                                
#                                 #get the dates into a list
#                                 date_in= open_f.ETo.lead.values
#                                 #get the start date
#                                 start_date = pd.to_datetime(open_f.S.values[0])
#                                 #Add the dates based on index value in date_in
#                                 date_out = []
#                                 for i in range(len(date_in)):
#                                     date_out.append(start_date + dt.timedelta(days = i))
                                
#                                 #Convert to julian date
#                                 end_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
                                
#                                 b_julian_out = [end_julian + i for i in range(len(date_out))]
                                
#                                 '''Find out if that file has a leap year, subtract appropriately'''
#                                 if pd.to_datetime(file[-14:-4]).year % 4 == 0:
#                                     subtract = 366
#                                     b_julian_out2 = [i-subtract if i>366 else i for i in b_julian_out]
#                                 else:
#                                     subtract = 365
#                                     b_julian_out2 = [i-subtract if i>365 else i for i in b_julian_out]
                                
                            
#                                 Et_ref_open_f = open_f.assign_coords(lead = b_julian_out2)
                                
#                                 '''Now we need to append to the dictionary with the same julian date values'''
#                                 for val in b_julian_out2:
#                                     try:
#                                         summation_ETo_mod0[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[0], X=i_X, Y=i_Y).values)})   
#                                         summation_ETo_mod1[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[1], X=i_X, Y=i_Y).values)})   
#                                         summation_ETo_mod2[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[2], X=i_X, Y=i_Y).values)})   
#                                         summation_ETo_mod3[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[3], X=i_X, Y=i_Y).values)})   

#                                     except KeyError:
#                                         pass
                                
#                         return(summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3)
                    
#                     #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
#                     ETo_7_day_mod0,ETo_7_day_mod1,ETo_7_day_mod2,ETo_7_day_mod3 = dict1_subx2()
                                        
#                     '''Now we have created a dictionary that contains the:
#                         1.) index value is the julian day
#                         2.) list of dictionaries containing:
#                         -init date file: summed 7-day value
#                         [{'1999-01-10': 20.289343},  {'1999-01-15': 25.726818}, ....}]

#                     Now the dictionary is filled with all values from all files, we need to sort
#                     each julian day by having the largest values as number 1 ranking. Then we can append the 
#                     Et_out file with the proper julian day and the EDDI value
                    
#                     EDDI calculation source :
#                         M. Hobbins, A. Wood, D. McEvoy, J. Huntington, C. Morton, M. Anderson, and C. Hain (June 2016): The Evaporative Demand Drought Index: Part I – Linking Drought Evolution to Variations in Evaporative Demand. J. Hydrometeor., 17(6),1745-1761, doi:10.1175/JHM-D-15-0121.1.
#                     '''
#                     def compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod):
#                         out_eddi_dictionary = {}
#                         #Key value in dictionary is julian_date
#                         for idx,julian_date in enumerate(ETo_7_day_mod):
                            
#                             out_eddi_dictionary[f'{julian_date}']=[]
                            
#                             subset_by_date = ETo_7_day_mod[f'{julian_date}']
                            
#                             #When looking at each julian date, now create a ranked list
#                             list_rank = []
#                             for date_init in range(len(subset_by_date)):
#                                 list_rank.append(list(subset_by_date[date_init].values())[0])
                            
#                             ''' i = 1 for maximum ETo aggregation in time window'''
#                             ranking = rankdata([-1 * i for i in list_rank]).astype(int)
#                             tukey = (ranking - 0.33)/ (len(ranking) + 0.33)
                                
#                             out_eddi = []
#                             '''Now we can calculate EDDI'''
#                             for probability in tukey:
#                                 if probability <= 0.5:
#                                     out_eddi.append(np.sqrt(-2*np.log(probability)))
#                                 else:
#                                     out_eddi.append((1-probability))  
                                        
#                             #constants
#                             c0 = 2.515517
#                             c1 = 0.802853
#                             c2 = 0.010328
#                             d1 = 1.432788
#                             d2 = 0.189269
#                             d3 = 0.001308
                                
#                             #Must reverse the sign of EDDI,so multiply by -1
#                             final_out_eddi = []
#                             for idx,w in enumerate(out_eddi):
#                                 final_out_eddi.append((w - ((c0 + c1*w + c2*w**2)/(1 + d1*w + d2*w**2 + d3*w**3)))*-1)
 
#                             out_eddi_dictionary[f'{julian_date}'].append(final_out_eddi)
                            
#                         return(out_eddi_dictionary)
                        
#                     EDDI_dict_mod0 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod = ETo_7_day_mod0)
#                     EDDI_dict_mod1 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod = ETo_7_day_mod1)
#                     EDDI_dict_mod2 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod = ETo_7_day_mod2)
#                     EDDI_dict_mod3 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod = ETo_7_day_mod3)

#                     '''Instead of re-looping through all of the files (very slow), we can start appending to 
#                     files one by one with the data that we have already collected.
                    
#                     The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
#                     has the initialized dates for each EDDI julian day summation'''
                    
#                     def improve_EDDI_dictionary(ETo_7_day_mod, EDDI_dict_mod):
#                         final_out_dictionary_all_eddi = {}

#                         for idx,julian_dattt in enumerate(ETo_7_day_mod):
#                             final_out_dictionary_all_eddi[f'{julian_dattt}'] = [] #create an empty list to append to
#                             sub_list = ETo_7_day_mod[f'{julian_dattt}'] #choose only the summation ETo with correct julian date
#                             sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)
    
#                             #sub list contains a dictionary for each julian date in current loop and the values of ETo
#                             #Save the init date values for each julian date
#                             for idxxx, init_date in enumerate(sub_list):
                                
#                                 sub_keys.append({list(init_date.keys())[0] :EDDI_dict_mod[f'{julian_dattt}'][0][idxxx]})
                            
                            
#                             final_out_dictionary_all_eddi[f'{julian_dattt}'].append(sub_keys) 
                            
#                         return(final_out_dictionary_all_eddi)
                    
#                     EDDI_next_dict_mod0 = improve_EDDI_dictionary(ETo_7_day_mod=ETo_7_day_mod0, EDDI_dict_mod=EDDI_dict_mod0)
#                     EDDI_next_dict_mod1 = improve_EDDI_dictionary(ETo_7_day_mod=ETo_7_day_mod1, EDDI_dict_mod=EDDI_dict_mod1)
#                     EDDI_next_dict_mod2 = improve_EDDI_dictionary(ETo_7_day_mod=ETo_7_day_mod2, EDDI_dict_mod=EDDI_dict_mod2)
#                     EDDI_next_dict_mod3 = improve_EDDI_dictionary(ETo_7_day_mod=ETo_7_day_mod3, EDDI_dict_mod=EDDI_dict_mod3)

                    
#                     '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
#                     we can append to aall EDDI files'''
                    
#                     #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
                
#                     '''Now that we have created new files, we can append each file with the data that was found'''
                    
#                     def add_to_npy_file(EDDI_next_dict,model_number):
                    
#                         for idx_,i_val in enumerate(EDDI_next_dict):
                          
#                         #for some reason it's a list in a list, this fixes that, loop through julian day
#                             EDDI_final_dict = EDDI_next_dict[f'{i_val}'][0]
                            
#                             for dic_init_and_eddi_val in EDDI_final_dict:
#                                 # print(dic_init_and_eddi_val)
#                                 #Open up the file and insert the value
#                                 init_day = list(dic_init_and_eddi_val.keys())[0]
                                
#                                 '''When appending using several scripts, sometimes the file will
#                                 load before another file has finished filling in the file and causing
#                                 a ValueError: cannot reshape array of size 15328 into shape (2,4,45,27,59)'''
 
#                                 try:
#                                     eddi_open = np.load(f'EDDI_{init_day}.npy',allow_pickle=True)
#                                 except ValueError:
#                                     os.system('sleep 2')
#                                     eddi_open = np.load(f'EDDI_{init_day}.npy',allow_pickle=True)
                                
#                                 lead_values = np.load(f'EDDI_{init_day}_julian_lead.npy',allow_pickle=True)
                                
#                                 #Sometimes we found indexes that really don't matter because they had to be 
#                                 try:
#                                     index_val=np.where(lead_values == int(i_val))[0][0]
#                                     eddi_open[0,model_number,index_val,i_Y,i_X] = list(dic_init_and_eddi_val.values())[0]
#                                     np.save(f'EDDI_{init_day}.npy',eddi_open)
#                                 except IndexError:
#                                     np.save(f'EDDI_{init_day}.npy',eddi_open)
#                                     pass
                        
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod0,model_number=0)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod1,model_number=1)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod2,model_number=2)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod3,model_number=3)
                        
                                                   
#                 #If not within CONUS, make np.nan
#                 else:
#                     #Open subX file for EDDI (initially it is empty with all 0s). If soil moisture SMERGE has no value, then EDDI has no value
                    
#                     try:
#                         eddi_file = np.load(f'EDDI_{_date}.npy')
#                     except ValueError:
#                         os.system('sleep 2')
#                         eddi_file[0,:,:,i_Y,i_X] = np.nan
#                         np.save(f'EDDI_{_date}.npy',eddi_file)      
 
#         print(f'Completed date {_date}')
          
# #%%

# if __name__ == '__main__':
#     p = Pool(num_processors)
#     p.map(multiProcess_EDDI_SubX_TEST, init_date_list[0:75])
   
###Old EDDI code, I just moved it into a new script for simplicity in looping
# #%%    
# '''process SubX files and create EDDI values'''
# # def multiProcess_EDDI_SubX_TEST(_date):
# for _date in init_date_list[0:12]:
#     try:
#         #don't re-work code if the netcdf file is already created
#         xr.open_dataset(f'{home_dir}/EDDI_{_date}.nc4')
#     except FileNotFoundError:
#         print(f'Calculating EDDI on SubX for {_date} with all models and saving as .npy in {home_dir}.')
# # for _date in init_date_list[0:80]:   
#         os.chdir(f'{home_dir}')
#         week_lead = 6

#         #Used for eliminating iterating over grid cells that don't matter
#         smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
      
    
#         #get julian day, timestamps, and the datetime
#         def date_file_info(_date, variable):
#             open_f = xr.open_dataset(f'{variable}_{_date}.nc4')
#             a_date_in= open_f[f'{variable}'].lead.values
#             #get the start date
#             a_start_date = pd.to_datetime(open_f.S.values[0])
#             #Add the dates based on index value in date_in
#             a_date_out = []
    
#             for a_i in range(len(a_date_in)):
#                 a_date_out.append(a_start_date + dt.timedelta(days = a_i))
                
#             start_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
#             #Julian day into a list                            
#             a_julian_out = [start_julian + i for i in range(len(a_date_out))]
            
#             #month of file
#             INdate_for_month = dt.datetime(int(_date[0:4]),int(_date[5:7]),int(_date[8:10]))
            
#             return(open_f,a_date_out,a_julian_out,INdate_for_month)
        
#         subx,file_timestamp_list,file_julian_list,file_datetime = date_file_info(_date=_date,variable='ETo')
         
    
#         #Now convert to julian date and append coordinates
#         subx2 = subx.assign_coords(lead = file_julian_list)

#         for i_Y in range(subx2.ETo.shape[3]):
#             for i_X in range(subx2.ETo.shape[4]):
#                 if _date == '1999-01-10':
#                     print(f'Working on lat {i_Y} and lon {i_X}')

                              
#                 '''Next if statement creates a mask to restrict within CONUS, because
#                 EDDI covers North America, this will create the mask using SMERGE plus 2 additional grid cell.
#                 Don't loop over the values that shouldn't have values.
#                 SMERGE indexes 38,6 and 38,7 are missing, but we can compute EDDI values '''
#                 if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

#                     def dict1_subx2():
#                         #append 7-day summation from all files to a new dictionary, 
#                         #trying to avoid having to loop through models.
#                         summation_ETo_mod0 = {}
#                         summation_ETo_mod1 = {}
#                         summation_ETo_mod2 = {}
#                         summation_ETo_mod3 = {}
#                         models = [0,1,2,3]
                        
#                         for julian_d in file_julian_list:
#                             #You must julian_d + week_lead because with EDDI you need 7-day looking backwards into the past. Must have 7 values.
#                             #Choose just model = 0 because we just need to know if there are 7-days total in any model
#                             if (len(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=0, X=i_X, Y=i_Y).values)) == 7:
#                                 summation_ETo_mod0[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod1[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod2[f'{week_lead+julian_d}']=[]
#                                 summation_ETo_mod3[f'{week_lead+julian_d}']=[]


#                                 summation_ETo_mod0[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[0], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod1[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[1], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod2[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[2], X=i_X, Y=i_Y).values)})   
#                                 summation_ETo_mod3[f'{week_lead+julian_d}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,week_lead+julian_d)).isel(S=0, model=models[3], X=i_X, Y=i_Y).values)})   

#                         '''7-day sum of Eto by index:
#                         Next we will append to each julian day value in the key in the dictionary with
#                         the same julian day from all files'''
                                
#                         for file in sorted(glob('ETo*.nc4')):
#                             #Remove unnecessary files that won't contain any of the dates, if difference in months
#                             #is >=2, then skip that file because of 45 day lead time.
                            
#                             OUTdate_for_month = dt.datetime(int(file[4:8]),int(file[9:11]),int(file[12:14]))
                                                        
#                             num_months = (OUTdate_for_month.month - file_datetime.month)
                            
#                             if num_months <3: #skip months that aren't within 3 months
#                                 #Dont' re-open the same file
#                                 if file[-14:-4] != _date:
                                    
#                                     open_f = xr.open_dataset(file)
                                    
#                                     '''Convert lead dates to a vector, then add it back into a netcdf
#                                     because I cannot convert np.datetime to a pd.datetime, I will need
#                                     to iterate over each of the dates in a list from Et_ref'''
                                    
#                                     #get the dates into a list
#                                     date_in= open_f.ETo.lead.values
#                                     #get the start date
#                                     start_date = pd.to_datetime(open_f.S.values[0])
#                                     #Add the dates based on index value in date_in
#                                     date_out = []
#                                     for i in range(len(date_in)):
#                                         date_out.append(start_date + dt.timedelta(days = i))
                                    
#                                     #Convert to julian date
#                                     end_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
                                    
#                                     b_julian_out = [end_julian + i for i in range(len(date_out))]
    
#                                     Et_ref_open_f = open_f.assign_coords(lead = b_julian_out)
                                    
#                                     '''Now we need to append to the dictionary with the same julian date values'''
#                                     for val in b_julian_out:
#                                         try:
#                                             summation_ETo_mod0[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[0], X=i_X, Y=i_Y).values)})   
#                                             summation_ETo_mod1[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[1], X=i_X, Y=i_Y).values)})   
#                                             summation_ETo_mod2[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[2], X=i_X, Y=i_Y).values)})   
#                                             summation_ETo_mod3[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=models[3], X=i_X, Y=i_Y).values)})   

#                                         except KeyError:
#                                             pass
                                    
#                         return(summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3)
                    
#                     #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
#                     ETo_7_day_mod0,ETo_7_day_mod1,ETo_7_day_mod2,ETo_7_day_mod3 = dict1_subx2()
                                        
#                     '''Now we have created a dictionary that contains the:
#                         1.) index value is the julian day
#                         2.) list of dictionaries containing:
#                         -init date file: summed 7-day value
#                         [{'1999-01-10': 20.289343},  {'1999-01-15': 25.726818}, ....}]

#                     Now the dictionary is filled with all values from all files, we need to sort
#                     each julian day by having the largest values as number 1 ranking. Then we can append the 
#                     Et_out file with the proper julian day and the EDDI value'''
                         
#                     def compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day):
#                         out_eddi_dictionary = {}
#                         for idx,julian_date in enumerate(ETo_7_day):
                            
#                             out_eddi_dictionary[f'{julian_date}']=[]
                            
#                             subset_by_date = ETo_7_day[f'{julian_date}']
                            
#                             #When looking at each julian date, now create a ranked list
#                             list_rank = []
#                             for date_init in range(len(subset_by_date)):
#                                 list_rank.append(list(subset_by_date[date_init].values())[0])
                            
                            
#                                 ranking = rankdata([-1 * i for i in list_rank]).astype(int)
#                                 tukey = (ranking - 0.33)/ (len(ranking) + 0.33)
                                
#                             out_eddi = []
#                             '''Now we can calculate EDDI'''
#                             for valu in tukey:
#                                 if valu <= 0.5:
#                                     w = out_eddi.append(np.sqrt(-2*np.log(valu)))
#                                 else:
#                                     w = out_eddi.append(1-valu)  
                                        
#                             #constants
#                             c0 = 2.515517
#                             c1 = 0.802853
#                             c2 = 0.010328
#                             d1 = 1.432788
#                             d2 = 0.189269
#                             d3 = 0.001308
                                
#                             final_out_eddi = []
#                             for w in out_eddi:
#                                 final_out_eddi.append(w - ((c0 + c1*w + c2*w**2)/(1 + d1*w + d2*w**2 + d3*w**3)))
                             
#                             out_eddi_dictionary[f'{julian_date}'].append(final_out_eddi)
                            
#                         return(out_eddi_dictionary)
                        
#                     EDDI_dict_mod0 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day = ETo_7_day_mod0)
#                     EDDI_dict_mod1 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day = ETo_7_day_mod1)
#                     EDDI_dict_mod2 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day = ETo_7_day_mod2)
#                     EDDI_dict_mod3 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day = ETo_7_day_mod3)

#                     '''Instead of re-looping through all of the files (very slow), we can start appending to 
#                     files one by one with the data that we have already collected.
                    
#                     The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
#                     has the initialized dates for each EDDI julian day summation'''
                    
#                     def improve_EDDI_dictionary(ETo_7_day_sum_dict, EDDI_dict):
#                         final_out_dictionary_all_eddi = {}

#                         for idx,julian_dattt in enumerate(ETo_7_day_sum_dict):
#                             final_out_dictionary_all_eddi[f'{julian_dattt}'] = [] #create an empty list to append to
#                             sub_list = ETo_7_day_sum_dict[f'{julian_dattt}'] #choose only the summation ETo with correct julian date
#                             sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)
    
#                             #sub list contains a dictionary for each julian date in current loop and the values of ETo
#                             #Save the init date values for each julian date
#                             for idxxx, dictonary in enumerate(sub_list):
                                
#                                 sub_keys.append({list(dictonary.keys())[0] :EDDI_dict[f'{julian_dattt}'][0][idxxx]})
                            
                            
#                             final_out_dictionary_all_eddi[f'{julian_dattt}'].append(sub_keys) 
                            
#                         return(final_out_dictionary_all_eddi)
                    
#                     EDDI_next_dict_mod0 = improve_EDDI_dictionary(ETo_7_day_sum_dict=ETo_7_day_mod0, EDDI_dict=EDDI_dict_mod0)
#                     EDDI_next_dict_mod1 = improve_EDDI_dictionary(ETo_7_day_sum_dict=ETo_7_day_mod1, EDDI_dict=EDDI_dict_mod1)
#                     EDDI_next_dict_mod2 = improve_EDDI_dictionary(ETo_7_day_sum_dict=ETo_7_day_mod2, EDDI_dict=EDDI_dict_mod2)
#                     EDDI_next_dict_mod3 = improve_EDDI_dictionary(ETo_7_day_sum_dict=ETo_7_day_mod3, EDDI_dict=EDDI_dict_mod3)

                    
#                     '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
#                     we can append to aall EDDI files'''
                    
#                     #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
                
#                     '''Now that we have created new files, we can append each file with the data that was found'''
                    
#                     def add_to_npy_file(EDDI_next_dict,model_number):
                    
#                         for idx_,i_val in enumerate(EDDI_next_dict):
                          
#                         #for some reason it's a list in a list, this fixes that, loop through julian day
#                             EDDI_final_dict = EDDI_next_dict[f'{i_val}'][0]
                            
#                             for dic_init_and_eddi_val in EDDI_final_dict:
#                                 # print(dic_init_and_eddi_val)
#                                 #Open up the file and insert the value
#                                 init_day = list(dic_init_and_eddi_val.keys())[0]
#                                 eddi_open = np.load(f'EDDI_{init_day}.npy',allow_pickle=True)
#                                 lead_values = np.load(f'EDDI_{init_day}_julian_lead.npy',allow_pickle=True)
    
#                                 try:
#                                     index_val=np.where(lead_values == int(i_val))[0][0]
#                                     eddi_open[0,model_number,index_val,i_Y,i_X] = list(dic_init_and_eddi_val.values())[0]
#                                     np.save(f'EDDI_{init_day}.npy',eddi_open)
#                                 except IndexError:
#                                     np.save(f'EDDI_{init_day}.npy',eddi_open)
#                                     pass
                        
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod0,model_number=0)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod1,model_number=1)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod2,model_number=2)
#                     add_to_npy_file(EDDI_next_dict = EDDI_next_dict_mod3,model_number=3)
                        
                                                   
#                 #If not within CONUS, make np.nan
#                 else:
#                     #Open subX file for EDDI (initially it is empty with all 0s). If soil moisture SMERGE has no value, then EDDI has no value
#                     eddi_file = np.load(f'EDDI_{_date}.npy')
#                     eddi_file[0,:,:,i_Y,i_X] = np.nan
#                     np.save(f'EDDI_{_date}.npy',eddi_file)      
 
#         print(f'Completed date {_date}')
          
# #%%
                   




# '''Convert each .npy file from multiProcess_EDDI_SubX output as a netcdf file'''
# #Read a subx file
# subX_file = xr.open_dataset('tasmax_GMAO_2015-09-12.nc')

# # test_load_a = f'{home_dir}/test_/EDDI_2005-05-25.npy'


# for date_1 in init_date_list:
    
#     test_load = np.load(f'{home_dir}/EDDI_{date_1}.npy')
#     # test_load = np.load(test_load_a)
    
    
#     #find S values
#     st_day =  pd.to_datetime(date_1)
#     st_day_2 = st_day + dt.timedelta(days=1)
#      # day_doyey = st_day.timetuple().tm_yday #julian day
    
#      #add more julian dates
#     day_julian_a = [pd.to_datetime(st_day_2) + dt.timedelta(days=i) for i in range(45)]
#     day_julian_b = [i.timetuple().tm_yday for i in day_julian_a]                            
    
#     S_values = [pd.to_datetime(date_1)+ dt.timedelta(days=1),pd.to_datetime(date_1)]
    
    
#     #convert to netcdf file
#     #Convert to an xarray object
#     var_OUT = xr.Dataset(
#         data_vars = dict(
#             EDDI = (['S','model','lead', 'Y','X'], test_load[:,:,:,:,:]),
#         ),
#         coords = dict(
#             S = S_values,
#             model = subX_file.M.values,
#             lead = day_julian_b,
#             Y = subX_file.Y.values,
#             X= subX_file.X.values,
#         ),
#         attrs = dict(
#             Description = 'Evaporative Demand Drought Index from SubX'),
#     )                    
    
#     #Save as a netcdf for later processing
#     var_OUT.to_netcdf(path = f'{home_dir}/EDDI_{date_1}.nc4', mode ='w')
#     print(f'Saved file into {home_dir}.')

#%%

'''This old code was a different formulation for reference ET. It didn't produce
meaningful results across CONUS (despite needing data more data as input than the current). 
This is the ASCE reference ET package'''
#  -- don't run    

# soil_moisture = xr.open_dataset(glob('mrso*')[0]).stack(grid = ['X','Y'])
# #Find nan locations for later use
# nan_dataset = np.where(np.isnan(soil_moisture.READ_DESCRIPTION_ATTRS.values), np.nan, 0)

# #Open file and change all values of 0.0 to nan
# temperature_1 = (xr.open_dataset(glob('tas*')[0]).stack(grid = ['X','Y']))
# temperature_1 = temperature_1.where(temperature_1 != 0)

# windV = (xr.open_dataset(glob('vas*')[0]).stack(grid = ['X','Y']))
# windV = windV.where(windV != 0)

# windU = (xr.open_dataset(glob('uas*')[0]).stack(grid = ['X','Y']))
# windU = windU.where(windU != 0)

# precip = (xr.open_dataset(glob('pr*')[0]).stack(grid = ['X','Y']))
# precip = precip.where(precip != 0)

# dewpoint_temp_1 = (xr.open_dataset(glob('tdps*')[0]).stack(grid = ['X','Y']))
# dewpoint_temp_1= dewpoint_temp_1.where(dewpoint_temp_1 != 0)

# shortwave_radiation = (xr.open_dataset(glob('dswrf*')[0]).stack(grid = ['X','Y']))
# shortwave_radiation = shortwave_radiation.where(shortwave_radiation != 0)

# elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False).stack(grid = ['lon','lat'])
# elevation = elevation.data[0,:]

# dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values
# #%%Temp is good, esat is good, delta looks good, 
# #Step 1, Compute ET0 using the lagged average ensemble

# #Constants
# Cn = 1600 #mm/d
# Cd = 0.38 #mm/s

# def convert_temperature(temperature, dewpoint_temp) -> float:
#     '''Convert temperature to Celsius (from Kelvin)'''
#     temperature = np.subtract(temperature_1, 273.15)
#     dewpoint_temp = np.subtract(dewpoint_temp_1, 273.15)
    
#     return(temperature, dewpoint_temp)

# temperature, dewpoint_temp = convert_temperature(temperature = temperature_1, \
#                                                  dewpoint_temp = dewpoint_temp_1)

# # #Check outputs
# # print(temperature.READ_DESCRIPTION_ATTRS.values)
# # print(temperature.READ_DESCRIPTION_ATTRS[0,200,0].values)
# # print(np.nanmax(temperature.READ_DESCRIPTION_ATTRS.values))
# # print(np.nanmin(temperature.READ_DESCRIPTION_ATTRS.values))

# # print(dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values)
# # print(dewpoint_temp.READ_DESCRIPTION_ATTRS[0,200,0].values)
# # print(np.nanmax(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
# # print(np.nanmin(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
# ##############################################################################
# test_temp = 20
# def esat_vapor_pressure(temperature, test_temp) -> float:
#     '''Tetens Formula: Need temp in Celsius. Returns values in kPa
#     Test temp should return 2339Pa or 2.34kPa'''
#     esat = 0.6112 * np.exp((17.27*temperature)/(237.3 + temperature))
#     esat_test = 0.6112 * np.exp((17.27*test_temp)/(237.3 + test_temp))
#     delta = (2504 * esat) / (temperature + 237.3)**2
#     return(esat, delta, esat_test)

# esat, delta, esat_test = esat_vapor_pressure(temperature = temperature, test_temp=test_temp)

# print()
# print(f'Esat_test kPa is {esat_test} at {test_temp} C.')
# print()
# # print(esat.READ_DESCRIPTION_ATTRS.values)
# print(f'Max esat value is {np.nanmax(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min esat value is {np.nanmin(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(delta.READ_DESCRIPTION_ATTRS.values)
# print(f'Max delta value is {np.nanmax(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min delta value is {np.nanmin(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# ##############################################################################
# test_dewpoint_temp = 16
# def e_vapor_pressure(esat, dewpoint_temp, test_dewpoint_temp) -> float:
#     '''Need temp in Celsius. Returns vapor pressure in Pa and humidity
#     values in percent (e.g., 0.78)
#     ea_test should return 1.81kPa'''
#     ea = 0.6108 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
#     ea_test = 0.6108 * np.exp((17.27*test_dewpoint_temp)/(237.3 + test_dewpoint_temp))
#     rel_humidity = np.divide(ea,esat)
#     return(ea, rel_humidity, ea_test)

# ea, rel_humidity, ea_test = e_vapor_pressure(esat = esat, dewpoint_temp = dewpoint_temp,test_dewpoint_temp = 16)

# print(f'Ea_test kPa is {ea_test} at {test_dewpoint_temp} C.')
# print()
# # print(ea.READ_DESCRIPTION_ATTRS.values)
# print(f'Max ea value is {np.nanmax(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min ea value is {np.nanmin(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(rel_humidity.READ_DESCRIPTION_ATTRS.values)
# print(f'Max rel_humidity value is {np.nanmax(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
# print()
# print(f'Min rel_humidity value is {np.nanmin(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
# print()
# ##############################################################################
# def convert_dswrf_RN(shortwave_radiation) -> float:
#     '''Returns shorwave radiation in MJ/m2 from W/m2'''
#     shortwave_rad = shortwave_radiation / 1000000
#     return(shortwave_rad)

# shortwave_rad = convert_dswrf_RN(shortwave_radiation = shortwave_radiation)

# # print(shortwave_rad.READ_DESCRIPTION_ATTRS.values)
# print(f'Max shortwave_rad value is {np.nanmax(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
# print()
# print(f'Min shortwave_rad value is {np.nanmin(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
# print()
# ##############################################################################
# #Square both u and v and then take the square root
# def compute_windSpeed(windU, windV) -> float:
#     ''' Returns average daily wind at 10-m height m/s'''
#     #Square u and v and take the square root to get windspeed
#     windSpeed = xr.ufuncs.sqrt(np.square(windU) + np.square(windV))
#     return (windSpeed)

# windSpeed = compute_windSpeed(windU = windU, windV = windV)

# # print(windSpeed.READ_DESCRIPTION_ATTRS.values)
# print(f'Max windSpeed value is {np.nanmax(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
# print()
# print(f'Min windSpeed value is {np.nanmin(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
# print()


# ##############################################################################
# def mean_atmospheric_pressure(elevation) -> float:
#     pressure = xr.zeros_like(temperature)

#     for k in range(elevation.data.shape[0]):
#         if k == 0 and i == 0 and j == 0:
#             print(f'Pressure at elevation {elevation[k].values} meters is {101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)} kPa.')
#         pressure.READ_DESCRIPTION_ATTRS[:,:,k] = 101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)
#         # pressure = pressure.where(temperature_1 != 0)

#     psychro = 0.000665 * pressure
#     '''Returns atmospheric pressure in kPa. And returns psychrometric constant'''
    
#     return(pressure, psychro)

# pressure, psychro = mean_atmospheric_pressure(elevation = elevation)


# # print(pressure.READ_DESCRIPTION_ATTRS.values)
# print()
# print(f'Max pressure value is {np.nanmax(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min pressure value is {np.nanmin(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(psychro.READ_DESCRIPTION_ATTRS.values)
# print(f'Max psychro value is {np.nanmax(psychro.READ_DESCRIPTION_ATTRS.values)}')
# print()
# print(f'Min psychro value is {np.nanmin(psychro.READ_DESCRIPTION_ATTRS.values)}')
# print()

# #%%
# #TODO
# '''NEED to convert all negative numbers to 0'''
# def compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, \
#                      esat, ea):
    
#     # ETref = xr.zeros_like(windSpeed)
#     '''Returns Reference Evapotranspiration in mm/d'''
#     ETref = ((0.408 * delta * shortwave_radiation) + (psychro*(Cn/temperature+273) * windSpeed * \
#         (esat - ea))) / (delta + (psychro * (1 + Cd * windSpeed)))
        
#     return(ETref)
    
# ETref_func = compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, esat, ea)

# ETref_func.READ_DESCRIPTION_ATTRS.values

# # print(ETref_func.READ_DESCRIPTION_ATTRS.values)
# print()
# print(f'Max ETref_func value is {np.nanmax(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
# print()
# print(f'Min ETref_func value is {np.nanmin(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
# print()


# ETref_unstack = ETref_func.unstack('grid').transpose('realization', 'time_series_dates', 'Y', 'X')

# #Convert to an xarray object
# var_OUT = xr.Dataset(
#     data_vars = dict(
#         ETo = (['realization','time_series_dates','Y','X'], ETref_unstack.READ_DESCRIPTION_ATTRS.values),
#     ),
#     coords = dict(
#         X = ETref_unstack.X.values,
#         Y = ETref_unstack.Y.values,
#         time_series_dates = ETref_unstack.time_series_dates.values,
#         realization = ETref_unstack.realization.values,
#     ),
#     attrs = dict(
#         Description = 'Reference crop evapotranspiration (mm/day)'),
# )

# var_OUT.ETo.values

# #Save as a netcdf for later processing
# var_OUT.to_netcdf(path = f'{home_dir}/ETo_lagged_average.nc4', mode ='w')
# #%%
# a = ETref_func.READ_DESCRIPTION_ATTRS[0,200,0].values
# print(a)


# #%%


# #%%

# for j in range(test.READ_DESCRIPTION_ATTRS.shape[1]):
#     for i in range(test.READ_DESCRIPTION_ATTRS.shape[2]):
#         if test.READ_DESCRIPTION_ATTRS[0,j,i] < 0:
#             print (test.READ_DESCRIPTION_ATTRS[0,i,i].values)
#             break
#         print (test.READ_DESCRIPTION_ATTRS[0,i,i].values) 
        
        
    

#%%
'''Test changing subX file order of dimensions to see if I can use cdo operators
to make a grid text file to re-grid.'''

# #Open a random subX file. Get grid description though CDO
# aa = a.transpose('lead','S', 'model','Y','X')
# bb = aa.rename_dims({'lead':'time'})

# bb.drop('S')

# cc = bb.ETo[:,:,:,:,:]



# cc = bb.drop_dims('S')
# dd = cc.drop_dims('model')



# #Convert to an xarray object
# var_OUT = xr.Dataset(
#     data_vars = dict(
#         ETo = (['time','Y','X'], bb.ETo[:,0,0,:,:]),
#     ),
#     coords = dict(
#         X = bb.X.values,
#         Y = bb.Y.values,
#         time = bb.lead.values),
#     attrs = dict(
#         Description = 'Reference crop evapotranspiration (mm/day)'),
# )                    

# #Save as a netcdf for later processing
# var_OUT.to_netcdf(path = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/test_/test_change_dim.nc4', mode ='w')
# print(f'Saved file into {home_dir}.')
