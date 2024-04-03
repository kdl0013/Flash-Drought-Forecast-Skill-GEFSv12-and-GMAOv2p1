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
from pyeto import fao



# dir1 = 'main_dir'
# num_processors = int('procs')
# mod = 'EMC'

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors = 7
mod = 'EMC'
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

init_date_list = [i[-14:-4] for i in sorted(glob('tasmax*.nc4'))]
if len(init_date_list) == 0:
    init_date_list = [i[-14:-4] for i in sorted(glob('tas*.nc4'))]


#%% GMAO
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX_GMAO(_date):

# for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')
        #Open up each file
        
        try:
            tasmax = np.subtract(xr.open_dataset(f'tasmax_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            tasmin = np.subtract(xr.open_dataset(f'tasmin_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        convert_W_to_MJ = 0.0864
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            dswrf = np.multiply(xr.open_dataset(f'dswrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass


        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            dlwrf = np.multiply(xr.open_dataset(f'dlwrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
            #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            ulwrf = np.multiply(xr.open_dataset(f'ulwrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass


        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            uswrf = np.multiply(xr.open_dataset(f'uswrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        try:
            tdps = np.subtract(xr.open_dataset(f'tdps_{mod}_{_date}.nc4'),273.15)
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
            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)
            windSpeed.to_dataset(name='windspeed').to_netcdf(f'windspeed_{mod}_{_date}.nc4')
            
            #net radiation
            short_rad = np.subtract(dswrf.dswrf,uswrf.uswrf).rename('short_rad')
            long_rad = np.subtract(dlwrf.dlwrf,ulwrf.ulwrf).rename('long_rad')
            
            srad = np.add(short_rad,long_rad).rename('srad').to_dataset()
            
            tavg = ((tasmax.tasmax+tasmin.tasmin)/2).to_dataset(name='tavg')

            try:
                srad = srad.assign_coords(S=np.atleast_1d(_date))
                srad.to_netcdf(f'srad_{mod}_{_date}.nc4')
            except ValueError:
                pass
            
            # tavg = (0.5*(tasmax.tasmax+tasmin.tasmin)).rename('tavg')*units.degC
            # tdps = tdps * units.degC
            
            # tdps.apply(fao.avp_from_tdew)
            # tdps.tdps.apply_ufunc( fao.avp_from_tdew)
            # # ea = 0.6112 * np.exp((17.27*tdps.tdps)/(237.3 + tdps.tdps))
            # fao.avp_from_tdew(tdps.tdps[0,0,0,10,10])
            
            # rel_humidity = xr.zeros_like(tdps.tdps)
            # for model_n in range(tasmax.M.shape[0]):
            #     for lead in range(tasmax.L.shape[0]):
            #         rel_humidity[0,model_n,lead,:,:] = metpy.calc.relative_humidity_from_dewpoint(temperature=tavg[0,model_n,lead,:,:], dewpoint=tdps.tdps[0,model_n,lead,:,:])
                
            # #Make mixing ratio
            # mix_ratio = xr.zeros_like(tdps.tdps)
            # for model_n in range(tasmax.M.shape[0]):
            #     for lead in range(tasmax.L.shape[0]):
            #         mix_ratio[0,model_n,lead,:,:] = metpy.calc.mixing_ratio_from_relative_humidity(pressure=atm_pressure.data[0,:,:], temperature=tavg[0,model_n,lead,:,:], relative_humidity=rel_humidity[0,model_n,lead,:,:])
                
            # #vapor pressure
            # ea = xr.zeros_like(tdps.tdps)
            # for model_n in range(tasmax.M.shape[0]):
            #     for lead in range(tasmax.L.shape[0]):
            #         ea[0,model_n,lead,:,:] = metpy.calc.vapor_pressure(pressure=atm_pressure.data[0,:,:], mixing_ratio=mix_ratio[0,model_n,lead,:,:])
    

            # ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')

            ea = xr.zeros_like(tdps.tdps)
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
            # fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values)
         
            for model_n in range(tasmax.M.shape[0]):
                for i_lead in range(tasmax.L.shape[0]):
                    for i_Y in range(tasmax.Y.shape[0]):
                        for i_X in range(tasmax.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                # print(i_X)
                                date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                                
                            # fao56_penman_monteith(net_rad=MJ/m2/d, t=Kelvin, ws=ms/s, svp=kPa, avp=kPa, delta_svp, psy, shf=0.0)
                                output_f.tasmin[0,model_n, i_lead, i_Y, i_X] = \
                                        fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                                                  t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                                                  ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                                                  svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values)))
                                            
                                            
                                # print(fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                #                           t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                #                           ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                #                           svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values))))
                                    
                                ea[0,model_n, i_lead, i_Y, i_X] = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values)
                                        
            ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')

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
            
            
#%% EMC
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX_EMC(_date):
    # _date=init_date_list[0]
# for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')
        #Open up each file

        try:
            tasmax = np.subtract(xr.open_dataset(f'tasmax_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            tasmin = np.subtract(xr.open_dataset(f'tasmin_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        
        convert_W_to_MJ = 0.0864
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            dswrf = np.multiply(xr.open_dataset(f'dswrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass


        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            dlwrf = np.multiply(xr.open_dataset(f'dlwrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
            #86000 seconds in 1 day, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            ulwrf = np.multiply(xr.open_dataset(f'ulwrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass


        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            uswrf = np.multiply(xr.open_dataset(f'uswrf_{mod}_{_date}.nc4'),convert_W_to_MJ)
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        try:
            # specific humidity
            huss = xr.open_dataset(f'huss_{mod}_{_date}.nc4')

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
            and ('huss' in list(locals().keys())) and ('windU' in list(locals().keys())) \
                and ('windV' in list(locals().keys())) and ('tasmax' in list(locals().keys())) \
                    and ('ulwrf' in list(locals().keys())) and ('uswrf' in list(locals().keys())) \
                        and ('dlwrf' in list(locals().keys()))):
            # print('calc')
            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)
            windSpeed.to_dataset(name='windspeed').to_netcdf(f'windspeed_{mod}_{_date}.nc4')

            #net radiation
            short_rad = np.subtract(dswrf.dswrf,uswrf.uswrf).rename('short_rad')
            long_rad = np.subtract(dlwrf.dlwrf,ulwrf.ulwrf).rename('long_rad')
            
            srad = np.add(short_rad,long_rad).rename('srad').to_dataset()
            try:
                srad = srad.assign_coords(S=np.atleast_1d(_date))
                srad.to_netcdf(f'srad_{mod}_{_date}.nc4')
            except ValueError:
                pass
            
            tavg = ((tasmax.tasmax+tasmin.tasmin)/2).to_dataset(name='tavg')
  
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
            
            ea = xr.zeros_like(tavg.tavg)

            
            for model_n in range(tasmax.M.shape[0]):
                for i_lead in range(tasmax.L.shape[0]):
                    for i_Y in range(tasmax.Y.shape[0]):
                        for i_X in range(tasmax.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                # print(i_X)
                                date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                                
                                rh = qair2rh(qair=huss.spfh[0,model_n, i_lead, i_Y, i_X].values, \
                                             press=fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values), \
                                                 temp = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values)
                                        
                                '''Can calculate vapor pressure using specific humidity. Convert specific humidity into
                                relative humidity. Calculate vapor pressure with relative humidity, vapor pressure of tmin and tmax'''
                                # fao56_penman_monteith(net_rad=MJ/m2/d, t=Kelvin, ws=ms/s, svp=kPa, avp=kPa, delta_svp, psy, shf=0.0)
                                output_f.tasmin[0,model_n, i_lead, i_Y, i_X] = \
                                        fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                                                  t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                                                  ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                                                  svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  avp = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(tasmin.tasmin[0,model_n, i_lead, i_Y, i_X].values), \
                                                                                            svp_tmax=fao.svp_from_t(tasmax.tasmax[0,model_n, i_lead, i_Y, i_X].values), \
                                                                                            rh_mean = rh),
                                                                  delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values)))
                                            
                                            
                                # print(fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                #                           t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                #                           ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                #                           svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values))))
                                    
                                ea[0,model_n, i_lead, i_Y, i_X] = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(tasmin.tasmin[0,model_n, i_lead, i_Y, i_X].values), \
                                                          svp_tmax=fao.svp_from_t(tasmax.tasmax[0,model_n, i_lead, i_Y, i_X].values), \
                                                          rh_mean = rh)
                                        
            ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')

    
    
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

#%% RSMAS
def multiProcess_Refet_SubX_RSMAS(_date):
    # _date=init_date_list[0]
    
# for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')

        try:
            tasmax = np.subtract(xr.open_dataset(f'tasmax_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            tasmin = np.subtract(xr.open_dataset(f'tasmin_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass            
        
        convert_W_to_MJ = 0.0864
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            rad = np.multiply(xr.open_dataset(f'rad_{mod}_{_date}.nc4'),convert_W_to_MJ)

        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        try:
            # specific humidity
            huss = xr.open_dataset(f'huss_{mod}_{_date}.nc4')

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
        if ('tasmin' in list(locals().keys())) and ('rad' in list(locals().keys())) \
            and ('huss' in list(locals().keys())) and ('windU' in list(locals().keys())) \
                and ('windV' in list(locals().keys())) and ('tasmax' in list(locals().keys())):
                    
            # print('calc')
            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)
            windSpeed.to_dataset(name='windspeed').to_netcdf(f'windspeed_{mod}_{_date}.nc4')

            srad = rad.rename(rad='srad')
            try:
                srad = srad.assign_coords(S=np.atleast_1d(_date))
                srad.to_netcdf(f'srad_{mod}_{_date}.nc4')
            except ValueError:
                pass

            tavg = (0.5*(tasmax.tasmax+tasmin.tasmin)).rename('tavg')
            
            def qair2rh (qair, temp, press):
                press = np.multiply(press,10) #convert from kPa to millibars
                es =  6.112 * np.exp((17.67 * temp)/(temp + 243.5))
                e = qair * press / (0.378 * qair + 0.622)
                rh = e / es
                rh[rh > 1] <- 1
                rh[rh < 0] <- 0
                return(rh)
            
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
            ea = xr.zeros_like(tavg)
       
            
            for model_n in range(tasmax.M.shape[0]):
                for i_lead in range(tasmax.L.shape[0]):
                    for i_Y in range(tasmax.Y.shape[0]):
                        for i_X in range(tasmax.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                # print(i_X)
                                date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                                
                                rh = qair2rh(qair=huss.huss[0,model_n, i_lead, i_Y, i_X].values, \
                                             press=fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values), \
                                                 temp = tavg[0,model_n, i_lead, i_Y, i_X].values)
                                        
                                '''Can calculate vapor pressure using specific humidity. Convert specific humidity into
                                relative humidity. Calculate vapor pressure with relative humidity, vapor pressure of tmin and tmax'''
                                # fao56_penman_monteith(net_rad=MJ/m2/d, t=Kelvin, ws=ms/s, svp=kPa, avp=kPa, delta_svp, psy, shf=0.0)
                                output_f.tasmin[0,model_n, i_lead, i_Y, i_X] = \
                                        fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                                                  t = tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                                                  ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                                                  svp = fao.svp_from_t(tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  avp = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(tasmin.tasmin[0,model_n, i_lead, i_Y, i_X].values), \
                                                                                            svp_tmax=fao.svp_from_t(tasmax.tasmax[0,model_n, i_lead, i_Y, i_X].values), \
                                                                                            rh_mean = rh),
                                                                  delta_svp = fao.delta_svp(tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values)))
                                            
                                            
                                # print(fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                #                           t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                #                           ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                #                           svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values))))
                                    
                                ea[0,model_n, i_lead, i_Y, i_X] = fao.avp_from_rhmean(svp_tmin=fao.svp_from_t(tasmin.tasmin[0,model_n, i_lead, i_Y, i_X].values), \
                                                          svp_tmax=fao.svp_from_t(tasmax.tasmax[0,model_n, i_lead, i_Y, i_X].values), \
                                                          rh_mean = rh)
                                        
            ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')
       
       
       
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
#%% ESRL and ECCC
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX_ESRL_ECCC(_date):
    # _date=init_date_list[0]
# for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')
        try:
            tasmax = np.subtract(xr.open_dataset(f'tasmax_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        try:
            tasmin = np.subtract(xr.open_dataset(f'tasmin_{mod}_{_date}.nc4'),273.15)#convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass            
        
        convert_W_to_MJ = 0.0864
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            rad = np.multiply(xr.open_dataset(f'rad_{mod}_{_date}.nc4'),convert_W_to_MJ)

        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
 
       
        try:
            tdps = np.subtract(xr.open_dataset(f'tdps_{mod}_{_date}.nc4'),273.15)
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
        if (('tasmin' in list(locals().keys())) and ('tasmax' in list(locals().keys())) \
            and ('tdps' in list(locals().keys())) and ('windU' in list(locals().keys())) \
                and ('windV' in list(locals().keys())) and ('rad' in list(locals().keys()))):
            # print('calc')
            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)
            windSpeed.to_dataset(name='windspeed').to_netcdf(f'windspeed_{mod}_{_date}.nc4')
            
            srad = rad.rename(rad='srad')
            try:
                srad = srad.assign_coords(S=np.atleast_1d(_date))
                srad.to_netcdf(f'srad_{mod}_{_date}.nc4')
            except ValueError:
                pass
            
            tavg = (0.5*(tasmax.tasmax+tasmin.tasmin)).rename('tavg')
            
            ea = xr.zeros_like(tdps.tdps)
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
            # fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values)
         
            for model_n in range(tasmax.M.shape[0]):
                for i_lead in range(tasmax.L.shape[0]):
                    for i_Y in range(tasmax.Y.shape[0]):
                        for i_X in range(tasmax.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                # print(i_X)
                                date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                                
                            # fao56_penman_monteith(net_rad=MJ/m2/d, t=Kelvin, ws=ms/s, svp=kPa, avp=kPa, delta_svp, psy, shf=0.0)
                                output_f.tasmin[0,model_n, i_lead, i_Y, i_X] = \
                                        fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                                                  t = tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                                                  ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                                                  svp = fao.svp_from_t(tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  delta_svp = fao.delta_svp(tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values)))
                                            
                                            
                                # print(fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                #                           t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                #                           ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                #                           svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values))))
                                    
                                ea[0,model_n, i_lead, i_Y, i_X] = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values)
                                        
            ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')

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
            

#%% NRL
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX_NRL(_date):
    # _date=init_date_list[0]
# for _date in init_date_list:
    # print(_date)
    fileOUT_name=f'ETo_Penman_{mod}_{_date}.nc4'
    try:
        xr.open_dataset(f'{fileOUT_name}')
        #To delete and re-run
        # xr.open_dataset('this.nc')
        # print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except FileNotFoundError:
        print(f'Working on date {mod} {_date} to calculate Penman-Monteith ETo.')
        try:
            tavg = np.subtract(xr.open_dataset(f'tas_{mod}_{_date}.nc4'),273.15) #convert to celsius
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass

        convert_W_to_MJ = 0.0864
        try:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            rad = np.multiply(xr.open_dataset(f'rad_{mod}_{_date}.nc4'),convert_W_to_MJ)

        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
 
       
        try:
            tdps = np.subtract(xr.open_dataset(f'tdps_{mod}_{_date}.nc4'),273.15*2) #there is an error with this dataset
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
        if (('tavg' in list(locals().keys())) \
            and ('tdps' in list(locals().keys())) and ('windU' in list(locals().keys())) \
                and ('windV' in list(locals().keys())) and ('rad' in list(locals().keys()))):
            # print('calc')
            def compute_windSpeed(windU, windV) -> float:
                ''' Returns average daily wind at 10-m height m/s'''
                #Square u and v and take the square root to get windspeed
                windSpeed = np.sqrt((np.square(windU.uas) + np.square(windV.vas)))
                return (windSpeed)
        
            windSpeed = compute_windSpeed(windU = windU, windV = windV)
            windSpeed.to_dataset(name='windspeed').to_netcdf(f'windspeed_{mod}_{_date}.nc4')
            
            srad = rad.rename(rad='srad')
            try:
                srad = srad.assign_coords(S=np.atleast_1d(_date))
                srad.to_netcdf(f'srad_{mod}_{_date}.nc4')
            except ValueError:
                pass

            
            ea = xr.zeros_like(tdps.tdps)
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
            # fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values)
         
            for model_n in range(tavg.M.shape[0]):
                for i_lead in range(tavg.L.shape[0]):
                    for i_Y in range(tavg.Y.shape[0]):
                        for i_X in range(tavg.X.shape[0]):
                            
                            #Don't calculate extra data
                            if conus_mask.CONUS_mask[0,i_Y,i_X].values == 1:
                                # print(i_X)
                                date_val = pd.to_datetime(tavg.S.values[0]) + dt.timedelta(days=i_lead)
                                day_doy = date_val.timetuple().tm_yday #julian day
                                
                            # fao56_penman_monteith(net_rad=MJ/m2/d, t=Kelvin, ws=ms/s, svp=kPa, avp=kPa, delta_svp, psy, shf=0.0)
                                output_f.tas[0,model_n, i_lead, i_Y, i_X] = \
                                        fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                                                  t = tavg.tas[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                                                  ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                                                  svp = fao.svp_from_t(tavg.tas[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  delta_svp = fao.delta_svp(tavg.tas[0,model_n, i_lead, i_Y, i_X].values), \
                                                                  psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values)))
                                            
                                            
                                # print(fao.fao56_penman_monteith(net_rad = srad.srad[0,model_n, i_lead, i_Y, i_X].values,\
                                #                           t = tavg.tavg[0,model_n, i_lead, i_Y, i_X].values + 273.15, \
                                #                           ws= windSpeed[0,model_n, i_lead, i_Y, i_X].values, \
                                #                           svp = fao.svp_from_t(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           avp = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           delta_svp = fao.delta_svp(tavg.tavg[0,model_n, i_lead, i_Y, i_X].values), \
                                #                           psy = fao.psy_const(fao.atm_pressure(elevation_dir.data[0,i_Y, i_X].values))))
                                    
                                ea[0,model_n, i_lead, i_Y, i_X] = fao.avp_from_tdew(tdps.tdps[0,model_n, i_lead, i_Y, i_X].values)
                                        
            ea.to_dataset(name='vapor_pressure').to_netcdf(f'actual_vapor_pressure_{mod}_{_date}.nc4')

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
    if mod == 'GMAO':
        p.map(multiProcess_Refet_SubX_GMAO, init_date_list)
    elif mod == 'EMC':
        p.map(multiProcess_Refet_SubX_EMC, init_date_list)
    elif mod == 'RSMAS':
        p.map(multiProcess_Refet_SubX_RSMAS, init_date_list)        
    elif mod == 'ESRL' or mod == 'ECCC':
        p.map(multiProcess_Refet_SubX_ESRL_ECCC, init_date_list)               
    elif mod == 'NRL':
        p.map(multiProcess_Refet_SubX_NRL, init_date_list)               

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
