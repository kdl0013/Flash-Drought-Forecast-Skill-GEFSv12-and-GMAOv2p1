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

dir1 = 'main_dir'
model_NAM = 'model_name'

home_dir = f'{dir1}/Data/SubX/{model_name}'
os.chdir(home_dir)
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM = 'GMAO'

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
missing_file_list_anomaly= []
missing_file_list_original = []

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
        missing_file_list_anomaly.append(file)
        
    open_mrso = xr.open_dataset(f'mrso_{model_NAM}_{file[-14:-4]}.nc') 
    if np.count_nonzero(np.isnan(open_mrso.mrso[0,:,:,:,:].values)) == 286740:
        missing_file_list_original.append(file)


missing_file_list_anomaly
missing_file_list_original

#%%           




