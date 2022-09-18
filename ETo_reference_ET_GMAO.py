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

init_date_list = [i[-14:-4] for i in sorted(glob('tas*.nc4'))]

#%%
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX(_date):
    # _date=init_date_list[0]
# for _date in init_date_list:
    fileOUT_name=f'ETo_{mod}_{_date}.nc4'
    try:
        # xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        os.system(f'rm {home_dir}/ETo_{mod}_{_date}.nc4')

        xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
         
        print(f'Working on date {mod} {_date} to calculate Priestley-Taylor ETo.')
        #Open up each file
        
        '''GMAO, we are given up- and downwelling short- and longwave radiation.
        
        To get net radiation (complete budget)  Rn = down-short minus up-short  (+) down-long minus down-upwave
        '''
        
        
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
                # dswrf.dswrf.values
                #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
                mult_ = 0.0864
                dswrf = np.multiply(dswrf,mult_) #convert to MJ/m2
                # dswrf.dswrf.values
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass
    
    
            try:
                dlwrf = xr.open_dataset(f'dlwrf_{mod}_{_date}.nc4')
                mult_ = 0.0864
                dlwrf = np.multiply(dlwrf,mult_) #convert to MJ/m2
                # dlwrf.dlwrf.values
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass
            
            try:
                ulwrf = xr.open_dataset(f'ulwrf_{mod}_{_date}.nc4')
                mult_ = 0.0864
                ulwrf = np.multiply(ulwrf,mult_) #convert to MJ/m2
                # dlwrf.dlwrf.values
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass

            try:
                uswrf = xr.open_dataset(f'uswrf_{mod}_{_date}.nc4')
                mult_ = 0.0864
                uswrf = np.multiply(uswrf,mult_) #convert to MJ/m2
                # dlwrf.dlwrf.values
            except IndexError:
                pass
            except ValueError:
                pass
            except FileNotFoundError:
                pass

        #Check if object is already in memory, some models won't have it
        if ('dlwrf' in list(locals().keys())) and ('dswrf' in list(locals().keys())) \
            and ('ulwrf' in list(locals().keys())) and ('uswrf' in list(locals().keys())) \
                and ('tavg' in list(locals().keys())):
            
            short_rad = np.subtract(dswrf.dswrf,uswrf.uswrf).rename('short_rad')
            long_rad = np.subtract(dlwrf.dlwrf,ulwrf.ulwrf).rename('long_rad')
            
            srad = np.add(short_rad,long_rad).rename('srad').to_dataset()
      
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
            
            output_f = xr.zeros_like(tavg)
            
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
                        ETo = (['S', 'M','L','Y','X'],  output_f[list(output_f.keys())[0]].values),
                    ),
                    coords = dict(
                        X = output_f.X.values,
                        Y = output_f.Y.values,
                        L = julian_list,
                        M = output_f.M.values,
                        S = output_f.S.values
                    ),
                    attrs = dict(
                        Description = 'Reference crop evapotranspiration (mm/day). Priestley-Taylor formula'),
                )              
            except AttributeError:
                
                #Convert to an xarray object
                var_OUT = xr.Dataset(
                    data_vars = dict(
                        ETo = (['S', 'M','L','Y','X'],  output_f[list(output_f.keys())[0]].values),
                    ),
                    coords = dict(
                        X = output_f.X.values,
                        Y = output_f.Y.values,
                        L = julian_list,
                        M = output_f.model.values,
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
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)
