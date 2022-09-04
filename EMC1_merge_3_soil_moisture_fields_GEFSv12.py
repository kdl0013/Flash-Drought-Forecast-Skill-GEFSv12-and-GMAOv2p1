#!/usr/bin/env python3

'''
Combine the three soil moisture fields (total 0-100cm). Other models have about this
layer

'''
import os
import datetime as dt
import numpy as np
import xarray as xr
from glob import glob
import pandas as pd
import datetime as dt


dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'

# #for linux
save_out_dir_GEFS = f'{dir1}/Data/SubX/EMC/'
os.chdir(save_out_dir_GEFS)

#get all files of a specific soil moisutre field
all_files1 = sorted(glob('soilw1*.nc4'))
all_files2 = sorted(glob('soilw2*.nc4'))
all_files3 = sorted(glob('soilw3*.nc4'))
# all_files4 = sorted(glob('soilw4*.nc'))
#%%
for idx,i in enumerate(zip(all_files1,all_files2,all_files3)):
    open_f1 = xr.open_dataset(i[0])
    open_f2 = xr.open_dataset(i[1])
    open_f3 = xr.open_dataset(i[2])
    

    #Add three levels together
    open_f1[list(open_f1.keys())[0]][:,:,:,:,:]=\
    (open_f1[list(open_f1.keys())[0]].values + \
    open_f2[list(open_f2.keys())[0]].values + \
    open_f3[list(open_f3.keys())[0]].values)
        
    #files have both np.nan and 9.9999e+20
    
    open_f1=open_f1.rename(soilw1='RZSM')
    open_f1 = open_f1.where(open_f1['RZSM'] != 9999.) 
    file_date=pd.to_datetime(open_f1.S.values)[0]
    
    #Now change the file julian date for anomaly
    
    def julian_date(open_f1):
        #Return julian date for anomaly calculation
        a_date_in= len(open_f1.L.values)
        #get the start date

        a_start_date = pd.to_datetime(open_f1.S.values[0]) 

        a_date_out=[]
        for a_i in range(a_date_in):
            a_date_out.append((a_start_date + np.timedelta64(a_i,'D')).timetuple().tm_yday)

        return(a_date_out)
    

    #Add julian date
    julian_list = julian_date(open_f1)

    open_f1['L']= ("L",julian_list)
    open_f1.close()
    
    open_f1.to_netcdf(f'soilw_bgrnd_{file_date.year}-{file_date.month:02}-{file_date.day:02}.nc4')
    # ff='soilw_bgrnd_2022-06-29.nc4'
    # xr.open_dataset(ff)
# for i in vars_to_process:
#     os.system(f'mkdir -p {out_dir}/{i}')
#%%

#Open file, find out if it has 2 values. If it does, take the 2nd value

# for var in vars_to_process:
#     # var=vars_to_process[0]
#     print(f'Working on variable {var} to remove S dimension.')
#     all_files = sorted(glob(f'{var}*.nc'))
#     for file in all_files:
#         open_f = xr.open_dataset(file)
#         varname=list(open_f.keys())[0]
#         if len(open_f[f'{varname}'].S.values) == 2:
#             break
#             open_f = open_f.dropna(dim='S',how='all')
#             open_f.close()
#             len(np.unique(open_f[f'{varname}'].values))
#             open_f.to_netcdf(f'{file}4')
#             os.system(f'rm {file}')
#             os.system(f'mv {file}4 {file}')

