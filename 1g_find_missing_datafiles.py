#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
I noticed that some files were blank after ETo and RZSM anomaly processing e.g., 2004-06-24.

But ETo and RZSM anomaly for 2001-06-24 is fine.
@author: kdl
"""

import xarray as xr
import numpy as np
from glob import glob
import os

# dir1 = 'main_dir'
# model_NAM = 'model_name'

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM = 'GMAO'

home_dir = f'{dir1}/Data/SubX/{model_NAM}'
os.chdir(home_dir)


#Get date list
test_var='ETo'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{home_dir}/{test_var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list() 

#%%
'''Now we want to open up each ETo- and RZSM anomaly file (also EDDI), to see
which ones have missing data'''
var = 'RZSM'

sorted(glob(f'{var}_anomaly_*.nc4'))

#RZSM_anomaly
missing_anomaly_RZSM= [] #subx missing anomaly data files
missing_original_RZSM = [] #mrso missing data files
missing_originial_RZSM_m3_m3 = [] #missing converted mrso m3/m3 files

if var == 'EDDI':
    file_list = f'{var}_*.nc4'
else:
    file_list = f'{var}_anomaly_*.nc4'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    if np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values)) == 286740:
        missing_anomaly_RZSM.append(file)
        
    open_mrso = xr.open_dataset(f'mrso_{model_NAM}_{file[-14:-4]}.nc') 
    if np.count_nonzero(np.isnan(open_mrso.mrso[0,:,:,:,:].values)) == 286740:
        missing_original_RZSM.append(file)

    open_m3_m3 = xr.open_dataset(f'SM_converted_m3_m3/SM_SubX_m3_m3_{file[-14:-4]}.nc4') 
    if np.count_nonzero(np.isnan(open_m3_m3.SM_SubX_m3_m3_value[0,:,:,:,:].values)) == 286740:
        missing_originial_RZSM_m3_m3.append(file)

missing_anomaly_RZSM
missing_original_RZSM
missing_originial_RZSM_m3_m3


#%%
var = 'ETo'

missing_anomaly_ETo = []
missing_original_ETo = []

if var == 'EDDI':
    file_list = f'{var}_*.nc4'
else:
    file_list = f'{var}_anomaly_*.nc4'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    if np.count_nonzero(np.isnan(open_f.ETo_anom[0,:,:,:,:].values)) == 286740:
        missing_anomaly_ETo.append(file)
        
    open_ETo = xr.open_dataset(f'{var}_{file[-14:-4]}.nc4') 
    if np.count_nonzero(np.isnan(open_ETo.ETo[0,:,:,:,:].values)) == 286740:
        missing_original_ETo.append(file)


missing_anomaly_ETo
missing_original_ETo
#%%           
#No anomalies for EDDI
var = 'EDDI'

missing_original_EDDI = []

if var == 'EDDI':
    file_list = f'{var}_*.nc4'
else:
    file_list = f'{var}_anomaly_*.nc4'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    if np.count_nonzero(np.isnan(open_f.EDDI[0,:,:,:,:].values)) == 286740:
        missing_original_EDDI.append(file)
        

missing_original_EDDI


