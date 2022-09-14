#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    
    EMC dswrf and dlwrf units = W/m2 . FLUX, not net -- need to calculate
    GMAO - W/m2 . works fine . FLUX, not net -- need to calculate
    RSMAS - W/m2 -- net surface radiation https://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.RSMAS/.CCSM4/.hindcast/
    ESRL - W/m2 -- net surface radiation
    
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
import pandas as pd
from glob import glob
import refet
from multiprocessing import Pool
from datetime import timedelta
import datetime as dt

# dir1 = 'main_dir'
# num_processors = int('procs')
# mod = 'model_name'

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = 10
mod = 'GMAO'
var = 'ETo'

home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#Additional datasets
elevation_dir = xr.open_dataset(f'{dir1}/Data/elevation/elev_regrid.nc',decode_times=False)
ptC_dir = f'{dir1}/Data/Priestley_Taylor_Makkink_evap_coeff_maps'
albedo_dir = xr.open_dataset(f'{dir1}/Data/MERRA2/albedo_netshortFlux_netdownFlux_merged.nc4')

###Files
file_list = os.listdir()

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of _datethe pre-processing that is 
#completed)
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'tas')    

#%%
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX(_date):
    # _date=init_date_list[0]
# for _date in init_date_list:
    fileOUT_name=f'ETo_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')

        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
         
        print(f'Working on date {mod} {_date} to calculate Priestley-Taylor ETo.')
        #Open up each file
        
        
        
#%%        
        if mod == "GMAO" or mod == "EMC":
            #These two models need to have net radiation computed
            #Missing data
            
            try:
                tavg = xr.open_dataset(f'tas_{mod}_{_date}.nc4')
                #convert to celsius
                tavg = np.subtract(tavg,273.15)
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass

                
            try:
                dswrf = xr.open_dataset(f'dswrf_{mod}_{_date}.nc4')
                dswrf.dswrf.values
                #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
                divisor = 0.0864
                dswrf = np.multiply(dswrf,divisor) #convert to MJ/m2
                dswrf.dswrf.values
                #Must make dswrf into mm/day by mulitplying by .408 https://www.fao.org/3/x0490e/x0490e08.htm
                # dswrf = np.multiply(dswrf,0.408)#convert to mm/d
                
                #Now multiply shortwave radiation by 1 - albedo to get net shortwave
                day_select = pd.to_datetime(list(dswrf.S.values)[0]).strftime("%Y-%m-%d")
                #We now have the date of the subx file  
                merra_days = list(pd.to_datetime(list(albedo_dir.ALBEDO.time.values)).strftime("%Y-%m-%d"))
                
                merra_days.index(day_select)
                
                dswrf = np.multiply(dswrf,(1 - albedo_dir.ALBEDO.isel(time=0).values))
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass
    
    
            try:
                dlwrf = xr.open_dataset(f'dlwrf_{mod}_{_date}.nc4')
                dlwrf.dlwrf.values
                divisor = 0.0864
                dlwrf = np.multiply(dlwrf,divisor) #convert to MJ/m2
                dlwrf.dlwrf.values
                # dlwrf = np.multiply(dlwrf,0.408)#convert to mm/d
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass
            
        elif mod == 'RSMAS' or mod == 'ESRL':
            srad = xr.open_dataset(f'rad_{mod}_{_date}.nc4')


        #Check if object is already in memory, some models won't have it
        if ('dlwrf' in list(locals().keys())) and ('dswrf' in list(locals().keys())):
            srad = np.add(dlwrf.dlwrf,dswrf.dswrf).to_dataset()
        elif 'dswrf' in (list(locals().keys())):
            srad = dswrf
        elif 'srad' in (list(locals().keys())):
            srad = srad


        ptC = xr.open_dataset(f'{ptC_dir}/pt_coeff_FINAL.nc')

        #convert to lead to julian day for later processing anomalies
        def date_file_info(SubX_file):
            try:
                a_date_in= SubX_file.L.values
            except AttributeError:
                a_date_in= SubX_file.lead.values
 
            #get the start date
            a_start_date = pd.to_datetime(SubX_file.S.values[0])
            a_date_out=[]
            for a_i in range(len(a_date_in)):
                a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)

            return(a_date_out)
        
        #Only process if we have temperature 
        if ('tavg' in list(locals().keys())) and ('srad' in list(locals().keys())):
            julian_list = date_file_info(tavg)
            
            output_f = xr.zeros_like(tavg)
            
            def priestley_taylor(tavg1,srad1,elevation1,ptC1):
                atm_pressure = 101.3 * ((293 - 0.0065 * elevation1)/293)**5.26 # kPa
                gamma = 0.000665 * atm_pressure # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
                delta = 4098 * (0.6108 * np.exp(17.27 * tavg1 / (tavg1  + 237.3))) / (tavg1  + 237.3)**2 # Slope saturated vapor pressure
                soil_heat_flux = np.zeros_like(delta)
                PET = ptC1*delta/(delta+gamma)*(srad1-soil_heat_flux)
                #check PET, add more parenthesis...may not 
                return PET
    
            for i_mod in range(tavg[list(tavg.keys())[0]].shape[1]):
                for i_lead in range(tavg[list(tavg.keys())[0]].shape[2]):
                    
                    #Just put the except blocks because of formatting between model files
                    try:
                        tavg1=tavg[list(tavg.keys())[0]].isel(L=i_lead,M=i_mod).values[0]
                    except ValueError:
                        tavg1=tavg[list(tavg.keys())[0]].isel(lead=i_lead,model=i_mod).values[0]

                    try:
                        srad1=srad[list(srad.keys())[0]].isel(L=i_lead,M=i_mod).values[0]
                    except ValueError:
                        srad1=srad[list(srad.keys())[0]].isel(lead=i_lead,model=i_mod).values[0]

                    elevation1=elevation_dir[list(elevation_dir.keys())[0]].isel(time=0).values
                    ptC1=ptC[list(ptC.keys())[0]].values
                    
                    # output_f[list(output_f.keys())[0]][0,i_mod, i_lead, :, :] = priestley_taylor(tavg1,srad1,elevation1,ptC1)*1000
                    output_f[list(output_f.keys())[0]][0,i_mod, i_lead, :, :] = priestley_taylor(tavg1,srad1,elevation1,ptC1)
    
            #Sometime infinity values mess up the saving of the file
            output_f = output_f.where(output_f.apply(np.isfinite)).fillna(0.0)
            
            try:
                #Convert to an xarray object
                var_OUT = xr.Dataset(
                    data_vars = dict(
                        ETo = (['S', 'model','lead','Y','X'],  output_f[list(output_f.keys())[0]].values),
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
            except AttributeError:
                
                #Convert to an xarray object
                var_OUT = xr.Dataset(
                    data_vars = dict(
                        ETo = (['S', 'model','lead','Y','X'],  output_f[list(output_f.keys())[0]].values),
                    ),
                    coords = dict(
                        X = output_f.X.values,
                        Y = output_f.Y.values,
                        lead = julian_list,
                        model = output_f.model.values,
                        S = output_f.S.values
                    ),
                    attrs = dict(
                        Description = 'Reference crop evapotranspiration (mm/day). Priestley-Taylor formula'),
                )              

                
            #Save as a netcdf for later processing
            var_OUT.to_netcdf(path = f'{home_dir}/ETo_{mod}_{_date}.nc', mode ='w')
            #compress
            os.system(f'ncks -4 -L 1 {home_dir}/ETo_{mod}_{_date}.nc {home_dir}/{fileOUT_name}')
            os.system(f'rm {home_dir}/ETo_{mod}_{_date}.nc')
            print(f'Saved {_date} into {home_dir}.')


#%%
'''Compute Reference ET for MERRA-2 data.

NOTES:
    Net shortwave radiation already computed. Don't need albedo to net radiation

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
