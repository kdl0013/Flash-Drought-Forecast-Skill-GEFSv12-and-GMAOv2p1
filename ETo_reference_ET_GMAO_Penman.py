#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    
Penman monteith formula for ETo. This is a test.
    


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
from metpy.units import units
import metpy.calc

# dir1 = 'main_dir'
# num_processors = int('procs')
# mod = 'model_name'

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = 5
mod = 'GMAO'
var = 'ETo'

home_dir = f'{dir1}/Data/SubX/{mod}'
script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

#Additional datasets
elevation_dir = xr.open_dataset(f'{dir1}/Data/elevation/elev_regrid.nc',decode_times=False)
atm_pressure = (101.3 * ((293 - 0.0065 * elevation_dir)/293)**5.26)*units.kPa # kPa --- LOOKS FINE AFTER INSPECTION

#CONUS mask
conus_mask = xr.open_dataset(f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4')

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of _datethe pre-processing that is 
#completed)

init_date_list = [i[-14:-4] for i in sorted(glob('tas_GMAO*.nc4'))]

#%%
'''Compute Reference ET for SubX data'''
# def multiProcess_Refet_SubX(_date):
    # _date=init_date_list[0]
for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        # print('alreadYc')
        #To delete and re-run
        # os.system(f'rm {home_dir}/ETo_{mod}_{_date}.nc4')

        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')
        #Open up each file
        
        '''GMAO, we are given up- and downwelling short- and longwave radiation.
        
        To get net radiation (complete budget)  Rn = down-short minus up-short  (+) down-long minus down-upwave
        '''
    

        try:
            tasmax = xr.open_dataset(f'tasmax_{mod}_{_date}.nc4')
            #convert to celsius
            tasmax = np.subtract(tasmax,273.15)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            tasmin = xr.open_dataset(f'tasmin_{mod}_{_date}.nc4')
            #convert to celsius
            tasmin = np.subtract(tasmin,273.15)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
            
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
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
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            dlwrf = xr.open_dataset(f'dlwrf_{mod}_{_date}.nc4')
            # dswrf.dswrf.values
            #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
            mult_ = 0.0864
            dlwrf = np.multiply(dlwrf,mult_) #convert to MJ/m2
            # dswrf.dswrf.values
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            ulwrf = xr.open_dataset(f'ulwrf_{mod}_{_date}.nc4')
            # dswrf.dswrf.values
            #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
            mult_ = 0.0864
            ulwrf = np.multiply(ulwrf,mult_) #convert to MJ/m2
            # dswrf.dswrf.values
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass


        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            uswrf = xr.open_dataset(f'uswrf_{mod}_{_date}.nc4')
            # dswrf.dswrf.values
            #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
            mult_ = 0.0864
            uswrf = np.multiply(uswrf,mult_) #convert to MJ/m2
            # dswrf.dswrf.values
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        try:
            tdps = xr.open_dataset(f'tdps_{mod}_{_date}.nc4')
            #convert to celsius
            tdps = np.subtract(tdps,273.15)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
   

        try:
            windU = xr.open_dataset(f'uas_{mod}_{_date}.nc4')
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        try:
            windV = xr.open_dataset(f'vas_{mod}_{_date}.nc4')
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
    
        #Need all of them in memory to proceed with calculation

        #Check if object is already in memory, some models won't have it
        if (('dswrf' in list(locals().keys())) and ('tasmin' in list(locals().keys())) \
            and ('tdps' in list(locals().keys())) and ('windU' in list(locals().keys())) \
                and ('windV' in list(locals().keys())) and ('tasmax' in list(locals().keys())) \
                    and ('ulwrf' in list(locals().keys())) and ('uswrf' in list(locals().keys())) \
                        and ('dlwrf' in list(locals().keys()))):
            # print('calc')
            #net radiation
            short_rad = np.subtract(dswrf.dswrf,uswrf.uswrf).rename('short_rad')
            long_rad = np.subtract(dlwrf.dlwrf,ulwrf.ulwrf).rename('long_rad')
            
            srad = np.add(short_rad,long_rad).rename('srad').to_dataset()
            srad = srad.assign_coords(S=np.atleast_1d(_date))
            
            tavg = (tasmax.tasmax-tasmin.tasmin).rename('tavg')*units.degC
            tdps = tdps * units.degC
            
            rel_humidity = xr.zeros_like(tdps.tdps)
            for model_n in range(tasmax.M.shape[0]):
                for Y in range(tasmax.Y.shape[0]):
                    for X in range(tasmax.X.shape[0]):

                        rel_humidity[0,model_n,:,Y,X] = metpy.calc.relative_humidity_from_dewpoint(temperature=tavg[0,model_n,:,Y,X], dewpoint=tdps.tdps[0,model_n,:,Y,X])
            
            #Make mixing ratio
            mix_ratio = xr.zeros_like(tdps.tdps)
            for model_n in range(tasmax.M.shape[0]):
                for Y in range(tasmax.Y.shape[0]):
                    for X in range(tasmax.X.shape[0]):

                        mix_ratio[0,model_n,:,Y,X] = metpy.calc.mixing_ratio_from_relative_humidity(pressure=atm_pressure.data[0,Y,X], temperature=tavg[0,model_n,:,Y,X], relative_humidity=rel_humidity[0,model_n,:,Y,X])
            
            #vapor pressure
            ea = xr.zeros_like(tdps.tdps)
            for model_n in range(tasmax.M.shape[0]):
                for Y in range(tasmax.Y.shape[0]):
                    for X in range(tasmax.X.shape[0]):

                        ea[0,model_n,:,Y,X] = metpy.calc.vapor_pressure(pressure=atm_pressure.data[0,Y,X], mixing_ratio=mix_ratio[0,model_n,:,Y,X])


            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)


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
    

            julian_list = date_file_info(tasmin)
            
            output_f = xr.zeros_like(tasmin)
    
            for i_mod in range(tasmax.M.shape[0]):
                for i_lead in range(tasmax.L.shape[0]):
                    #changed these 2 lines of code from under i_X
      
                    for i_Y in range(tasmax.Y.shape[0]):
                        for i_X in range(tasmax.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                    
                                # print(i_X)
                                date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                            
                                output_f.tasmin[0,i_mod, i_lead, i_Y, i_X] = \
                                    (refet.Daily(tmin = tasmin.tasmin[0,i_mod, i_lead, i_Y, i_X].values, \
                                                tmax = tasmax.tasmax[0,i_mod, i_lead, i_Y, i_X].values, \
                                                ea = ea[0,i_mod, i_lead, i_Y, i_X].values, \
                                                rs = srad.srad[0,i_mod, i_lead, i_Y, i_X].values, \
                                                uz = windSpeed[0,i_mod, i_lead, i_Y, i_X].values, \
                                                zw = 10, elev = elevation_dir.data[0,i_Y, i_X].values,
                                                lat = tasmax.Y[i_Y].values,
                                                doy = day_doy, method = 'asce',
                                                input_units={'tmin': 'C', 'tmax': 'C', \
                                                             'rs': 'Langleys', 'uz': 'm/s', \
                                                                 'lat': 'deg'}).etr())[0]   
                                        

    
    
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
                    Description = 'Reference crop evapotranspiration (mm/day). Penman-Monteith formula'),
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