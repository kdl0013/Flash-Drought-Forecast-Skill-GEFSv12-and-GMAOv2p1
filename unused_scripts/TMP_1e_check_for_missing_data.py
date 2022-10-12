#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
   Find any missing files from SubX model for ETo & RZSM anomalies, and EDDI

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob

import warnings
warnings.filterwarnings("ignore")

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'

home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
script_dir = f'{dir1}/Scripts'
anom_dir = f'{home_dir}/anomaly'

os.chdir(home_dir)

gridMET_dir = f'{dir1}/Data/gridMET'
smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask file
###Files
file_list = os.listdir()

global var
var='ETo'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{home_dir}/{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list()   


'''Now we want to open up each ETo- and RZSM anomaly file (also EDDI), to see
which ones have missing data'''
# var = 'RZSM'

for var in ['RZSM','ETo']:
    print(f'Checking variable {var} for missing anomaly data in files.')
    file_list = f'{anom_dir}/{var}_anomaly_*.nc4'
    
    file_list_out = []
    for file in sorted(glob(file_list)):
        # print()
        if "LEAD" not in file:
            file_list_out.append(file)

        
            
    missing_anomaly_var = [] #subx missing anomaly data files
    missing_original_var = []
    
    # test_file = f'{home_dir}/test/ETo_anomaly_2000-09-17.nc'
    
    for file in file_list_out:
        open_f = xr.open_dataset(file)
        # open_f = xr.open_dataset(test_file)
    
        open_f.close()
        '''after further inspection, if a file has all nan values then it has a total
        of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
        therefore, this is a missing data file.
        '''
        
        #anomaly. I inspected a file that was filled correctly and it had 528 missing value
        #If file is completely empty, the file will have all zeros and a unique value of only 1
        
        for lead_values in range(open_f[f'{var}_anom'].shape[2]):
            if len(np.unique(open_f[f'{var}_anom'][:,:,lead_values,:,:])) <= 1:
                missing_anomaly_var.append(file)    
    
        '''I manually calculated ETo from SubX variables, look if I am missing any
        data from those files'''
    
        if var == 'ETo':
            open_ETo = xr.open_dataset(f'{var}_{file[-14:-4]}.nc') 
            if np.count_nonzero(np.isnan(open_ETo.ETo[0,:,:,:,:].values)) == 286740:
                missing_original_var.append(file)
        elif var == 'RZSM':
            open_m3_m3 = xr.open_dataset(f'SM_converted_m3_m3/SM_SubX_m3_m3_{file[-14:-4]}.nc4') 
            if np.count_nonzero(np.isnan(open_m3_m3.SM_SubX_m3_m3_value[0,:,:,:,:].values)) == 286740:
                missing_original_var.append(file)
    
    
    # missing_anomaly_var[0][-14:-4]
    missing_dates = [i[-14:-4] for i in missing_anomaly_var]
    
    print(f'Missing data in {len(missing_dates)}.')
    print(f'Missing {var} SubX data (no anomaly) in {len(missing_original_var)} files.')
