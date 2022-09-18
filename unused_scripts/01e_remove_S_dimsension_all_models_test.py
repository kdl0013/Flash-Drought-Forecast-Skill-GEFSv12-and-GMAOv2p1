#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove unneeded S dimension in SubX files.
Don't run multiprocessing, for some reason it was just messing up files. 
Couldn't figure out why

Also, switch coordinates using nco operators to be on same grid as EMC GEFSv12 and
MERRA 2.

Also add the actual julian day of the year for the file to assist with anomaly
calculation.

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob
import pandas as pd
import datetime as dt
from datetime import timedelta
from multiprocessing import Pool


# dir1 = 'main_dir'
# mod = 'model_name'
num_processors=8
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='GMAO'
subX_dir = f'{dir1}/Data/SubX/fromCasper'
out_dir = f'{subX_dir}/S_dim_removed'
os.system(f'mkdir -p {out_dir}')
os.chdir(subX_dir)

#%%#SMERGE Soil Moisture
#SubX.dat depth file

# if mod == 'RSMAS' or mod == 'ESRL':
#     variables  = ['pr','rad','tasmax','tasmin','mrso']
# elif mod == 'GMAO':
#     variables  = ['pr','dswrf','tasmax','tasmin','mrso']
# elif mod == 'EMC':
#     variables  = ['pr','soilw1','soilw2','soilw3','soilw4','tasmax','tasmin','dswrf'] 

variables  = ['pr','soilw1','soilw2','soilw3','soilw4','tas','dswrf','dlwrf','rad','mrso'] 

#test GMAO var dlwrf, not writing correctly
# var='dlwrf'
for var in variables:
    print(f'Working on variable {var}.')
# def remove_S_dim(var):
    
    sub1_dir = subX_dir
    var_name=var
    os.chdir(sub1_dir)
    # try:
    #     xr.open_dataset(sorted(glob(f'{var}*.nc4'))[0])
    # except IndexError:
   
    for file in sorted(glob(f'{var}*.nc')):
        try:
            # xr.open_dataset(f'{out_dir}/{file}4')
            xr.open_dataset(f'{out_dir}/making_new_files.nc')
        except FileNotFoundError:
            open_f = xr.open_dataset(file)
            out_test1 = np.empty_like(open_f.to_array()).squeeze()
            out_test2 = np.empty((1,out_test1.shape[1],out_test1.shape[2],\
                                out_test1.shape[3],out_test1.shape[4])) #Keep this dataset, 1 dimension
            '''RSMAS - Both S dimensions have values with RSMAS
                       S[0] dim = initialized date (name of file) + 7 days
                       S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                       
               EMC  - S[0] dim = initialized date (name of file) + 1 day - no data
                      S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                      
               GMAO  - S[0] dim = initialized date (name of file) + 1 day - no data
                      S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                      
               ESRL  - Both S dimensions have values with ESRL
                      S[0] dim = initialized date (name of file) + 7 day
                      S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                      
            '''
            
            # if ('GMAO' in file):
            #     out_test2 = out_test1[0,:,:,:,:]               
            # else:
            out_test2 = out_test1[1,:,:,:,:]
                
            
            def julian_date(open_f,file):
                #Return julian date for anomaly calculation
                a_date_in= len(open_f.L.values)
                #get the start date
                # if 'GMAO' in file:
                #     a_start_date = pd.to_datetime(open_f.S.values[0]) - np.timedelta64(1,'D')
                # else:
                a_start_date = pd.to_datetime(open_f.S.values[0]) 
    
                a_date_out=[]
                for a_i in range(a_date_in):
                    a_date_out.append((a_start_date + np.timedelta(days=a_i)).timetuple().tm_yday)
        
                return(a_date_out)
            

            #Add julian date
            julian_list = julian_date(open_f,file)
            out_test2
            open_f['L']= ("L",julian_list)
            
            #Drop S dimension to save storage space
            open_f = open_f.dropna(dim='S',how='all')
            
            # open_f = open_f.where(open_f.apply(np.isfinite)).fillna(0.0)
            # open_f=open_f.rename(S='time')
            
            open_f.close()
            # os.system(f'rm a_{file}4')
          
            #TODO: Change the lead times to be julian day of year for anomaly
            #FOR GMAO only, subtract 1 from the filename because it is 1 day off
    
            # open_ft=open_f.assign_coords({'test':julian_list})
            # open_ft.to_netcdf(f'{out_dir}/a_{file}')  #make a file, just do it

            
            open_f.to_netcdf(f'{out_dir}/a_{file}')  #make a file, just do it
            try:
                # julian_list = julian_date(open_f,file)
                
                if 'pr' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            pr = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'mrso' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            RZSM = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'dswrf' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            dswrf = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'dlwrf' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            dlwrf = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'rad' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            rad = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'tas' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            tas = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'soilw1' in file:
                    #Convert to an xarray object to save dates as 
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            soilw1 = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'soilw2' in file:
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            soilw2 = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'soilw3' in file:
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            soilw3 = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                elif 'soilw4' in file:
                    var_OUT = xr.Dataset(
                        data_vars = dict(
                            soilw4 = (['S', 'model','lead','Y','X'], open_f[list(open_f.keys())[0]].values),
                        ),
                        coords = dict(
                            X = open_f.X.values,
                            Y = open_f.Y.values,
                            L = julian_list,
                            M = open_f.M.values,
                            S = open_f.S.values
                        ),
                        attrs = dict(
                            Description = f'{open_f[list(open_f.keys())[0]].long_name} {open_f[list(open_f.keys())[0]].level_type} {open_f[list(open_f.keys())[0]].units}.'),
                    )   
                    
                    
                var_OUT.close()
                # var_OUT['L'] = ("L",julian_list) #doesn't work
                '''Won't save the julian date in the lead time...This is necessary
                for anomaly calculation'''
                
                var_OUT.to_netcdf(f'{out_dir}/a_{file}', mode='w', format='netcdf4', engine='netcdf4') #remake the file with new data
                
                #Flip the X,Y coordinates to match other datasets. 0,0 at top left corner
                if ('GMAO' in file) or ('ESRL' in file) or ('RSMAS' in file):
                    os.system(f'ncpdq -O -a -Y,X {out_dir}/a_{file} {out_dir}/b_{file}') #rename out_dir file
                    os.system(f'ncks -4 -L 1 {out_dir}/b_{file} {out_dir}/{file}4')
                    os.system(f'rm {out_dir}/a_{file} {out_dir}/b_{file}')
                else:
                    os.system(f'ncks -4 -L 1 {out_dir}/a_{file} {out_dir}/{file}4')
                    os.system(f'rm {out_dir}/a_{file}')
                        
            except IndexError:
                open_f.to_netcdf(f'{out_dir}/{file}4') #No data in files
            
            #re-assign coords (This doesn't work when saving for some dumb reason)
            # open_f=open_f.assign_coords({"L":julian_list})
            # open_f.L.values
            
            #Flip the X,Y coordinates to match other datasets. 0,0 at top left corner


#%%
# if __name__ == '__main__':
#     p=Pool(num_processors)
#     p.map(remove_S_dim,variables)

# variables = ['dswrf','tasmax', 'tasmin', 'uas', 'vas', 'mrso','cape','pr','tdps','SM_SubX']
# var = variables[0]

# #Add a new dataset for RZSM and save into home_dir
# sub1_dir = f'{subX_dir}/SM_converted_m3_m3'
# os.chdir(sub1_dir)

# for file in glob('*.nc4'):
#     try:
#         xr.open_dataset(f'{subX_dir}/RZ{file}')
#     except FileNotFoundError:
        
#         SubX_file=xr.open_dataset(file)
        
#         def date_file_info(SubX_file):
    
#             a_date_in= SubX_file.lead.values
#             #get the start date
#             a_start_date = pd.to_datetime(SubX_file.S.values[0])
#             a_date_out=[]
#             for a_i in range(len(a_date_in)):
#                 a_date_out.append((a_start_date + timedelta(days=a_i)).timetuple().tm_yday)
    
#             return(a_date_out)
    
#         julian_list = date_file_info(SubX_file)
        
#         #re-assign coords
#         SubX_file=SubX_file.assign_coords(lead=julian_list)
        
#         #save to home directory
#         SubX_file.to_netcdf(f'{subX_dir}/RZ{file}')
