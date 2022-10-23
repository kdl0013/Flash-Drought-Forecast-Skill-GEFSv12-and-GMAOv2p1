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

SEARCH
Convert numpytimedelta to regular timestamp
dDateTime = dt.datetime.utcfromtimestamp(a.tolist()/1e9)

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
num_processors=9
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod='model_name'
subX_dir = f'{dir1}/Data/SubX/fromCasper'
out_dir = f'{subX_dir}/S_dim_removed'
os.system(f'mkdir -p {out_dir}')
os.chdir(subX_dir)

# cp_output = f'{dir1}/Data/SubX/{mod}
# def return_output_to_new_folder():
#     os.system(c'')
#remove all files
# os.system(f'rm {out_dir}/*{mod}*')

# glob('*NRL*')
#%%#SMERGE Soil Moisture
#SubX.dat depth file

# if mod == 'RSMAS' or mod == 'ESRL':
#     variables  = ['pr','rad','tasmax','tasmin','mrso']
# elif mod == 'GMAO':
#     variables  = ['pr','dswrf','tasmax','tasmin','mrso']
# elif mod == 'EMC':
#     variables  = ['pr','soilw1','soilw2','soilw3','soilw4','tasmax','tasmin','dswrf'] 

variables  = ['huss','pr','soilw1','soilw2','soilw3','soilw4','tas','dswrf','dlwrf','rad','mrso','ulwrf','uswrf','tasmin','tasmax','tdps','uas','vas'] 
# variables  = ['huss'] 

#test GMAO var dlwrf, not writing correctly
# var='dlwrf'
# for var in variables:

def remove_S_dim(var):
    print(f'Working on variable {var}.')
    sub1_dir = subX_dir
    var_name=var
    os.chdir(sub1_dir)
    # try:
    #     xr.open_dataset(sorted(glob(f'{var}*.nc4'))[0])
    # except IndexError:
        
    def julian_date(open_f,file):
        '''Day of year'''
        #Return julian date (doy) for anomaly calculation
        a_date_in= len(open_f.L.values)
        #get the start date
        a_start_date =pd.to_datetime(file.split('_')[-1].split('.')[0])

        a_date_out=[]
        for a_i in range(a_date_in):
            a_date_out.append((a_start_date + np.timedelta64(a_i,'D')).timetuple().tm_yday)

        return(a_date_out)
    #%%
    def name_and_return_files(file,open_f):
        #Returns the new name of the file
        if 'pr' in file:
            #Convert to an xarray object to save dates as 
            var_OUT = xr.Dataset(
                data_vars = dict(
                    pr = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    RZSM = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    dswrf = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    dlwrf = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    rad = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    tas = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    soilw1 = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    soilw2 = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    soilw3 = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
                    soilw4 = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'ulwrf' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    ulwrf = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'uswrf' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    uswrf = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'tasmin' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    tasmin = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
            
        elif 'tasmax' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    tasmax = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'uas' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    uas = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'vas' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    vas = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
    
        elif 'huss' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    huss = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        elif 'tdps' in file:
            var_OUT = xr.Dataset(
                data_vars = dict(
                    tdps = (['S', 'M','L','Y','X'], open_f[list(open_f.keys())[0]].values),
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
        
        return(var_OUT)
    
    
#%%    
    #Remove old file
    os.system(f'rm {out_dir}/*.tmp')
    for file in sorted(glob(f'{var}*{mod}*.nc')):
        try:
            xr.open_dataset(f'{out_dir}/{file}4')
            # xr.open_dataset(f'no_file.nc')
        except FileNotFoundError:
        #%%
        # var='tdps'
        # sub1_dir = subX_dir
        # var_name=var
        # os.chdir(sub1_dir)
        # for file in sorted(glob(f'{var}*{mod}*.nc')):
            
            if mod == 'NRL':
                '''This model only has 1 realization. So just need to rename and then change some things for anomaly calculation'''
                open_f = xr.open_dataset(file)
                julian_list = julian_date(open_f,file)
                open_f['L']= ("L",julian_list)
                open_f['S'] = np.atleast_1d(file.split('.')[0].split('_')[-1])
                
                open_f=name_and_return_files(file, open_f)
                open_f.to_netcdf(f'a_{file}') #save temp file, need to reorient
                #Now reorient
                os.system(f'ncpdq -O -a -Y,X a_{file} {out_dir}/{file}4') #rename out_dir file
            else:
                open_f = xr.open_dataset(file)

                if var=='huss':
                    #need to remove the P layer
                    open_f=open_f.drop('P')
                    open_f = open_f.huss[0,:,:,:,:].to_dataset()
    
    
                # out_test = np.empty_like(open_f.to_array()).squeeze()
                '''RSMAS - Both S dimensions have values with RSMAS
                           S[0] dim = initialized date (name of file) + 7 days
                           S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                           
                   EMC  - S[0] dim = initialized date (name of file) + 1 day - no data
                          S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                          
                   GMAO  - S[0] dim = initialized date (name of file) + 1 day - no data
                          S[1] dim = initialized date (name of file) - KEEP THIS DIM
                          
                   ESRL  - Both S dimensions have values with ESRL
                          S[0] dim = initialized date (name of file) + 7 day
                          S[1] dim = initialized date (name of file) --- KEEP THIS DIM
                 
                   ECCC  - S[0] dim = day of file
                          S[1] dim = initialized date (with data) (name of file) --- KEEP THIS DIM
                 
                   NRL - S[0] dim is the only dim
                 
                    
                '''

                
                #Remove S dim
                open_f[f'{var_name}'][0,:,:,:,:] = np.nan
                #Add julian date
                julian_list = julian_date(open_f,file)
                open_f['L']= ("L",julian_list)
                #Drop S dimension (needed for future assistance with anomaly skill)
                open_f = open_f.dropna(dim='S',how='all')
                open_f = name_and_return_files(file,open_f)

                # open_f = open_f.where(open_f.apply(np.isfinite)).fillna(0.0)
                # open_f=open_f.rename(S='time')
                
                open_f.close()
                # os.system(f'rm a_{file}4')
    
                open_f.to_netcdf(f'{out_dir}/a_{file}')  #make a file, just do it
                
                try:
                    # julian_list = julian_date(open_f,file)
     
    
    
                    #Flip the X,Y coordinates to match other datasets. 0,0 at top left corner
                    if ('GMAO' in file) or ('ESRL' in file) or ('RSMAS' in file) or ('ECCC' in file):
                        os.system(f'ncpdq -O -a -Y,X {out_dir}/a_{file} {out_dir}/{file}4') #rename out_dir file
                        # os.system(f'ncks -4 -L 1 {out_dir}/b_{file} {out_dir}/{file}4')
                        os.system(f'rm {out_dir}/a_{file}')
                    else:
                        # os.system(f'ncks -4 -L 1 {out_dir}/a_{file} {out_dir}/{file}4')
                        os.system(f'mv {out_dir}/a_{file} {out_dir}/{file}4')
                            
                except IndexError:
                    pass
                    # open_f.to_netcdf(f'{out_dir}/{file}4') #No data in files
                
                #re-assign coords (This doesn't work when saving for some dumb reason)
                # open_f=open_f.assign_coords({"L":julian_list})
                # open_f.L.values
                
                #Flip the X,Y coordinates to match other datasets. 0,0 at top left corner



#%%
if __name__ == '__main__':
    p=Pool(num_processors)
    p.map(remove_S_dim,variables)

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
