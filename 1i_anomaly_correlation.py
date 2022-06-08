#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 19:15:03 2022

@author: kdl
"""

import xarray as xr
import numpy as np
import numpy.ma as ma
import os
import datetime as dt
import pandas as pd
from glob import glob
from scipy.stats import rankdata
import sys

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
mod = 'GMAO'
subX_dir = f'{dir1}/Data/SubX/{mod}/' #where subX model data lies 
mask_file = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4'


conus_mask = xr.open_dataset(mask_file) #open mask file
#High plains CONUS mask
HP_conus_mask = conus_mask['USDM-HP_mask']
#West CONUS_mask
West_conus_mask = conus_mask['USDM-West_mask']


#Test variable
var='EDDI'

os.chdir(subX_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob(f'{var}*_*-*.nc4'))
date_list = [i[-14:-4] for i in date_list]


#EDDI actual values
eddi_subX_dir = f'{dir1}/Data/EDDI/EDDI_SubX_values' 

#%%
'''Now find the anomaly correlation by lead time'''

eddi_subx_all = xr.open_mfdataset(f'{subX_dir}/{var}*.nc4')
eddi_obs_all = xr.open_mfdataset(f'{eddi_subX_dir}/{var}*.nc')

# eddi_subx_all.EDDI[0,:,6,15,15].values #lead week 1

start_lead = 6
lead=0
lead_weeks_as_index = np.arange(start_lead,42,7)

week_1 = eddi_subx_all.EDDI[0,:,start_lead+lead,:,:].where(HP_conus_mask == 1)

xr.corr(week_1.where(HP_conus_mask==1), week_1.where(HP_conus_mask==1)).mean().compute()

eddi_subx_all[0]

xr.corr(eddi_subx_all.where(HP_conus_mask == 1),eddi_subx_all.where(HP_conus_mask == 1)).mean().

eddi_ = eddi_subx_all.where(HP_conus_mask == 1).mean().compute()







subX_dir = f'{dir1}/Data/SubX/{mod}'

subX_SM_out_dir = f'{subX_dir}/SM_converted_m3_m3' #conversion
#gridMET
gridMET_dir = f'{dir1}/Data/gridMET' #reference evapotranspiration
output_dir = f'{gridMET_dir}/ETo_SubX_values' #Refernce ET output directory
#RZSM
smerge_in_dir = f'{dir1}/Data/SMERGE_SM/Raw_data' #raw files that have been boxed in CONUS
SM_SubX_out_dir = f'{dir1}/Data/SMERGE_SM/Raw_data/SM_SubX_values' #Smerge values overlayed on SubX grid
#EDDI
eddi_dir = f'{dir1}/Data/EDDI/convert_2_nc'
eddi_subX_dir = f'{dir1}/Data/EDDI/EDDI_SubX_values'

image_dir = f'{dir1}/Outputs/MSE_plots'




#%%



















#%%
'''Since I've already created the actual values from SubX files for each index 
(EDDI, ETo, and RZSM) [from script 1d_RZSM_and_scatterplots.py], I can 
just re-factor portions of that code to calculate the pearson's correlation coefficent'''


















def flatten_RZSM_ETo_EDDI(date_list):
    '''Open RZSM anomalies and SubX file'''
    
    #Flatten 1 file to see the length
    ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().values.ravel())    
    file_length = len(date_list)
    final_length = ravel_length * file_length
    
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_RZSM = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open SMERGE file and ravel, and append to empty array
        smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].values.ravel()
        #append smerge to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_RZSM[count_value:ravel_end,0] = smerge_file[:]
        correlation_RZSM[count_value:ravel_end,1] = sub_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
 
    corr_df = pd.DataFrame(correlation_RZSM)
    corr_df.corr()
    
    '''Open gridMET and SubX file'''
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_ETo = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open gridMET file and ravel, and append to empty array
        g_file = xr.open_dataset(f'{output_dir}/ETo_SubX_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        s_file = xr.open_dataset(f'{subX_dir}/ETo_{date_}.nc4').ETo[0,:,:,:,:].values.ravel()
        #append gridMET to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_ETo[count_value:ravel_end,0] = g_file[:]
        correlation_ETo[count_value:ravel_end,1] = s_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
        
    corr_df = pd.DataFrame(correlation_ETo)
    corr_df.corr()   
    
    '''Open EDDI and SubX file'''
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_EDDI = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open EDDI file and ravel, and append to empty array
        g_file = xr.open_dataset(f'{eddi_subX_dir}/EDDI_SubX_{mod}_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        s_file = xr.open_dataset(f'{subX_dir}/EDDI_{date_}.nc4').EDDI[0,:,:,:,:].values.ravel()
        #append gridMET to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_EDDI[count_value:ravel_end,0] = g_file[:]
        correlation_EDDI[count_value:ravel_end,1] = s_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
    
    '''Currently very low correlation'''
    corr_df = pd.DataFrame(correlation_EDDI)
    corr_df.corr()
          
    return(correlation_RZSM, correlation_ETo, correlation_EDDI)