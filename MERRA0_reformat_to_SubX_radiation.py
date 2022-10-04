#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script will find create a new dataset of actual observations, but with
the format of SubX data. This will make correlation between datasets easier 
to calculate.

FOR GMAO GEOS-5, file name is technically day 0.

"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from multiprocessing import Pool

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod='GMAO'
home_dir = f'{dir1}/Data/SubX/{mod}'
obs_dir = f'{dir1}/Data/MERRA2' #reference evapotranspiration

#Open merra file, create radiation
merra_file = xr.open_dataset(f'{obs_dir}/radiation.nc4')
radiation_merra = np.add(merra_file.LWGNET,merra_file.SWNETSRF).rename('srad').to_dataset()
radiation_merra.to_netcdf(f'{obs_dir}/radiation_total.nc4')
n_processors = 10


#CONUS mask
conus_mask = xr.open_dataset(f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4')

os.chdir(home_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob(f'tas_{mod}*.nc4'))
date_list = [i[-14:-4] for i in date_list]


# _date = date_list[0]
# _date='1999-12-06'
#%%
def save_radiation(_date) -> float:    
    var='radiation'

    obs_file_name = f'{var}_{mod}_{_date}.nc4'
    
    '''First get the actual radiation values and save into a file. Create from SubX'''

    try:
        xr.open_dataset(f'{home_dir}/{obs_file_name}')
        # xr.open_dataset('test.nc') #only to make new files
        print(f'Already completed date {_date}. Saved in {home_dir}')
    except FileNotFoundError:
        print(f'Working on initialized day {_date} to find MERRA values from SubX models, leads, & coordinates and saving data into {home_dir}.')

        
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
            #86000 seconds in 1 home_dirday, 1000000J in 1 MJ https://www.anycodings.com/1questions/2366349/solar-energy-conversion-wm2-to-mjm2
            mult_ = 0.0864
            uswrf = np.multiply(uswrf,mult_) #convert to MJ/m2
            # dswrf.dswrf.values
        except IndexError:
            pass
        except ValueError:
            pass
        except FileNotFoundError:
            pass
        
        if (('dswrf' in list(locals().keys())) and ('ulwrf' in list(locals().keys())) and \
            ('uswrf' in list(locals().keys())) and ('dlwrf' in list(locals().keys()))):
            
            short_rad = np.subtract(dswrf.dswrf,uswrf.uswrf).rename('short_rad')
            long_rad = np.subtract(dlwrf.dlwrf,ulwrf.ulwrf).rename('long_rad')
            
            srad = np.add(short_rad,long_rad).rename('srad').to_dataset()
            srad=srad.assign_coords(S=np.atleast_1d(pd.to_datetime(_date)))
            srad=srad.assign_coords(L=np.arange(srad.L.shape[0]))
            
            srad.to_netcdf(f'{home_dir}/{obs_file_name}')

    return(0)

#%%
def convert_radiation_OBS_to_SubX(_date) -> float:   
    var='radiation'
    date_list = sorted(glob(f'{var}_{mod}*.nc4'))
    date_list = [i[-14:-4] for i in date_list]
    

    output_obs_dir = f'{obs_dir}/{var}_SubX_values' #Refernce ET output directory
    os.system(f'mkdir -p {output_obs_dir}')
    obs_file_name = f'{var}_SubX_{mod}_{_date}.nc4'
    # os.system(f'mkdir -p {output_obs_dir}')
    obs_radiation = xr.open_dataset(f'{obs_dir}/radiation_total.nc4')
    
    try:
        xr.open_dataset(f'{output_obs_dir}/{obs_file_name}')
        # xr.open_dataset('test.nc') #only to make new files
        print(f'Already completed date {_date}. Saved in {output_obs_dir}')
    except FileNotFoundError:
        print(f'Working on initialized day {_date} to find MERRA values from SubX models, leads, & coordinates and saving data into {output_obs_dir}.')
        
        
        #Open up anomaly file to get proper formatting
        sub_file = xr.open_dataset(f'{home_dir}/{var}_{mod}_{_date}.nc4')
        out_file = xr.zeros_like(sub_file)
    
        #Open gridMET and subset to the correct dates as SubX for full time series
        '''Issue with the dates'''
        dates_list = list(obs_radiation.time.values)
        df = pd.DataFrame({'dates':dates_list})
        df['dates'] = pd.to_datetime(df['dates'])
        df['dates'] = df['dates'].dt.floor('d')
        
        new_date_list = np.array(df['dates'])
        
        obs_radiation=obs_radiation.assign_coords(time=new_date_list)
        
        output_f = xr.zeros_like(sub_file)
 
        for i_lead in range(sub_file[list(sub_file.keys())[0]].L.shape[0]):
            #Because the S value is actually the first day of the forecast, don't add one
            date_val = pd.to_datetime(pd.to_datetime(_date) + dt.timedelta(days=i_lead))
            #1 day appears to have not been calculated because of julian days and 
            #anomaly code for +/-42 day (specically December 31, 2000)
            #Just take the average of the other two days before and after
            if np.count_nonzero(obs_radiation[list(obs_radiation.keys())[0]].sel(time = date_val).values == 0) == 1593:
                date_val1 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead-1)
                date_val2 = pd.to_datetime(sub_file.S.values[0]) + dt.timedelta(days=i_lead+1)
                
                output_f[list(output_f.keys())[0]][0,:, i_lead, :, :] = \
                np.nanmean( obs_radiation[list(obs_radiation.keys())[0]].sel(time = date_val1).values + \
                     obs_radiation[list(obs_radiation.keys())[0]].sel(time = date_val2).values)
            else:
                output_f[list(output_f.keys())[0]][0,:, i_lead, :, :] = \
                    obs_radiation[list(obs_radiation.keys())[0]].sel(time = date_val).values
                
                
        #Only keep 1 model
        out_file = output_f.isel(M=0)
        #TODO: Only keep the first 7 leads (total of 6 weeks). 0 index is 12 hour lead from initialization.
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                radiation = (['S','L','Y','X'],   out_file[list(output_f.keys())[0]].values),
            ),
            coords = dict(
                S = np.atleast_1d(pd.to_datetime(_date)),
                X = sub_file.X.values,
                Y = sub_file.Y.values,
                L = np.arange(sub_file.L.shape[0]),
        
            ),
            attrs = dict(
                Description = f'MERRA2 {var} values on the exact same date and grid \
                cell as {mod} SubX data'),
        )                    
        
 
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{output_obs_dir}/{obs_file_name}', mode ='w')
        
                # '''Compress to save space since I do not have to write again to file'''
                # os.system(f'ncks -4 -L 1 {output_ETo_dir}/{file_name} {output_ETo_dir}/eto_test.nc4')
                # os.system(f'mv {output_ETo_dir}/eto_test.nc4 {output_ETo_dir}/{file_name}')
 
        print(f'Completed ETo_SubX_{_date}')  
        
     

    return(0)
#%%
if __name__ == '__main__':
    p = Pool(n_processors)
    p.map(save_radiation,date_list)
    p.map(convert_radiation_OBS_to_SubX,date_list)
#%%
