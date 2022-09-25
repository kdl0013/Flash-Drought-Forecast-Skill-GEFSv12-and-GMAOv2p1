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

    def priestley_taylor(tavg1,srad1,elevation1,ptC1):
        atm_pressure = 101.3 * ((293 - 0.0065 * elevation1)/293)**5.26 # kPa --- LOOKS FINE AFTER INSPECTION
        specific_heat = 0.001013 #MJ/kg C
        mle_ratio_weight_vapor = 0.622
        latent_heat_of_vaporization = 2.501 - (0.002361 * tavg1) #MJ/kg, needs temp in Celsius
        
        '''gamma practically same results as old_gamma variable which uses only 0.000665 as the constant'''
        gamma = np.divide(np.multiply(atm_pressure,specific_heat), (mle_ratio_weight_vapor*latent_heat_of_vaporization)) #kPA
        # old_gamma = np.multiply(0.000665, atm_pressure) # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
        
        delta = 4098 * (0.6108 * np.exp(17.27 * tavg1 / (tavg1  + 237.3))) / (tavg1  + 237.3)**2 # Slope saturated vapor pressure #kPA
        soil_heat_flux = np.zeros_like(delta)
        PET = (ptC1*(delta/(delta+gamma))*(srad1-soil_heat_flux))/latent_heat_of_vaporization
        
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
mod = 'RSMAS'
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
    fileOUT_name=f'ETo_Priestley_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')

        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
         
        print(f'Working on date {mod} {_date} to calculate Priestley-Taylor ETo.')
        #Open up each file
        
        
        
#%%        
        
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
            srad = xr.open_dataset(f'rad_{mod}_{_date}.nc4')#In w/m2 already
            # srad.rad[0,0,0,10,10].values
            srad = np.multiply(srad,0.0864) #convert to MJ/m2
            # srad.rad[0,0,0,10,10].values

        except IndexError:
            #If file doesn't exist
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass



        #Check if object is already in memory, some models won't have it. Can only
        #complete ETo if we have temperature and radiation
        if ('tavg' in list(locals().keys())) and ('srad' in list(locals().keys())):

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
            

            julian_list = date_file_info(tavg)
            
            save_output = xr.zeros_like(srad)

            
            def priestley_taylor(tavg1,srad1,elevation1,ptC1):
                atm_pressure = 101.3 * ((293 - 0.0065 * elevation1)/293)**5.26 # kPa --- LOOKS FINE AFTER INSPECTION
                specific_heat = 0.001013 #MJ/kg C
                mle_ratio_weight_vapor = 0.622
                latent_heat_of_vaporization = 2.501 - (0.002361 * tavg1) #MJ/kg, needs temp in Celsius
                
                '''gamma practically same results as old_gamma variable which uses only 0.000665 as the constant'''
                gamma = np.divide(np.multiply(atm_pressure,specific_heat), (mle_ratio_weight_vapor*latent_heat_of_vaporization)) #kPA
                # old_gamma = np.multiply(0.000665, atm_pressure) # Psychrometric constant Cp/(2.45 * 0.622 = 0.000665
                
                delta = 4098 * (0.6108 * np.exp(17.27 * tavg1 / (tavg1  + 237.3))) / (tavg1  + 237.3)**2 # Slope saturated vapor pressure #kPA
                soil_heat_flux = np.zeros_like(delta)
                PET = (ptC1*(delta/(delta+gamma))*(srad1-soil_heat_flux))/latent_heat_of_vaporization
    
                return PET

            for i_mod in range(tavg[list(tavg.keys())[0]].M.shape[0]):
                for i_lead in range(tavg[list(tavg.keys())[0]].L.shape[0]):
                    
                    tavg1=tavg[list(tavg.keys())[0]].isel(L=i_lead,M=i_mod).values[0]

                    srad1=srad[list(srad.keys())[0]].isel(L=i_lead,M=i_mod).values[0]

                    elevation1=elevation_dir[list(elevation_dir.keys())[0]].isel(time=0).values
                    ptC1=ptC[list(ptC.keys())[0]].values
                    
                    # output_f[list(output_f.keys())[0]][0,i_mod, i_lead, :, :] = priestley_taylor(tavg1,srad1,elevation1,ptC1)*1000
                    save_output[list(save_output.keys())[0]][0,i_mod, i_lead, :, :] = priestley_taylor(tavg1,srad1,elevation1,ptC1)


            
            try:
                #Convert to an xarray object
                var_OUT = xr.Dataset(
                    data_vars = dict(
                        ETo = (['S', 'M','L','Y','X'],  save_output[list(save_output.keys())[0]].values),
                    ),
                    coords = dict(
                        X = save_output.X.values,
                        Y = save_output.Y.values,
                        L = julian_list,
                        M = save_output.M.values,
                        S = save_output.S.values
                    ),
                    attrs = dict(
                        Description = 'Reference crop evapotranspiration (mm/day). Priestley-Taylor formula'),
                )

            except AttributeError:
                
                #Convert to an xarray object
                var_OUT = xr.Dataset(
                    data_vars = dict(
                        ETo = (['S', 'M','L','Y','X'],  save_output[list(save_output.keys())[0]].values),
                    ),
                    coords = dict(
                        X = save_output.X.values,
                        Y = save_output.Y.values,
                        L = julian_list,
                        M = save_output.model.values,
                        S = save_output.S.values
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
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)
