#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Soil moisture percentile drop (SMPD) = decrease from above the 40th percentile to below the 20th percentile in 
3 weeks or less.

-I am going to use RZSM anomalies because it has already integrated the weekly mean.
If I used regular mrso files, 

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob
from scipy.stats import percentileofscore as pos
from multiprocessing import Pool


# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'
# var = 'RZSM'
# num_processors = int(7)

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'GMAO'
num_processors = int('7')

subx_RZSM_dir = f'{dir1}/Data/SubX/{model_NAM1}/SM_converted_m3_m3'

anom_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
percentile_dir = f'{anom_dir}/percentiles_anom_RZSM' 

os.system(f'mkdir -p {percentile_dir}')
var = 'RZSM'

#%% Create datasets for anomaly percentiles
def get_files_and_data_setup(anom_dir,subx_RZSM_dir):
    os.chdir(f'{anom_dir}') 
    
    global file_list_anomaly
    file_list_anomaly = sorted(glob('RZSM_anomaly*.nc4'))
    # file=file_list_anomaly[0]
    
    all_anom = xr.open_mfdataset(file_list_anomaly, concat_dim=['S'], combine='nested')
    
    global all_anomaly_numpy
    all_anomaly_numpy = all_anom.RZSM_anom.to_numpy()
    
    global file_list_RZSM
    os.chdir(f'{subx_RZSM_dir}')
    file_list_RZSM = sorted(glob('SM_SubX_**.nc4'))
    
    all_rzsm = xr.open_mfdataset(file_list_RZSM, concat_dim=['S'], combine='nested')
    global all_rzsm_numpy 
    all_rzsm_numpy = all_rzsm.SM_SubX_m3_m3_value.to_numpy()
    
    return(file_list_anomaly,all_anomaly_numpy,file_list_RZSM,all_rzsm_numpy)

'''Important note:
    Because the anomalies are positive for lower values, we must subract the 
    output by 100 to be properly organized'''
def create_anomaly_percentiles(file):
    
    try:
        xr.open_dataset(f'{percentile_dir}/anomaly_percentiles_{file[-14:-4]}.nc4')
    except FileNotFoundError:
        print(f'Working on anomaly percentiles {file}')      
        os.chdir(f'{anom_dir}') 
        # file=file_list_anomaly[0]
        open_f = xr.open_dataset(file)
        
        sm_input = open_f.RZSM_anom.to_numpy()
        sm_output = np.zeros_like(open_f.RZSM_anom)
        sm_output = sm_output[:,:,0:7,:,:] #only want first 7 leads
        
        def find_percentile_of_score(sm_input, sm_output):
            for model in range(sm_input.shape[1]):
                for Y in range(sm_input.shape[3]):
                    for X in range(sm_input.shape[4]):
                        for lead in range(7): #only looking at first 6 weeks for GMAO
                            all_values =all_anomaly_numpy[:,model,lead,Y,X]
                            # all_values = get_full_distribution(model,lead,Y,X)
                            
                            all_values = sorted(all_values)
                            all_values = [i for i in all_values if i !=0.0]
                            
                            #get percentile ranking (subtract 100 because higher anomalies are located in the upper percentiles)
                            sm_output[0,model,lead,Y,X] = 100 -pos(all_values,sm_input[:,model,lead,Y,X])
                        #if rzsm percentile is below the 20th: perform code else week is not in flash drought
            return(sm_output)
        
        sm_output = find_percentile_of_score(sm_input, sm_output)
        
        #Testing
        # model=3
        # lead=6
        # Y=17
        # X=44
        sm_input[0,3,:,17,44]
        sm_output[0,3,:,17,44]
                
        sm_input[0,3,:,10,10]
        sm_output[0,3,:,10,10]
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                Anomaly_percentile = (['S','model','lead','Y','X'], sm_output[:,:,0:7,:,:]),
            ),
            coords = dict(
                X = open_f.X.values,
                Y = open_f.Y.values,
                lead = np.arange(0,7),
                S = open_f.S.values,
                model = open_f.model.values
            ),
            attrs = dict(
                Description = 'Soil Moisture RZSM percentiles from anomaly data.'),
        )                    
    

        var_final.to_netcdf(path = f'{percentile_dir}/anomaly_percentiles_{file[-14:-4]}.nc4')
        var_final.close()
        
        return()

    
#%% Create datasets for RZSM volumetric water content percentiles


'''Important note:
    Because the anomalies are positive for lower values, we must subract the 
    output by 100 to be properly organized'''
def create_RZSM_percentiles(file):
    
    try:
        xr.open_dataset(f'{percentile_dir}/rzsm_percentiles_{file[-14:-4]}.nc4')
    except FileNotFoundError:
        print(f'Working on rzsm percentiles {file}')      
        #must change to directory to open file due to structure of list
        os.chdir(f'{subx_RZSM_dir}') 
        # file=file_list_anomaly[0]
        open_f = xr.open_dataset(file)
        
        sm_input = open_f.SM_SubX_m3_m3_value.to_numpy()
        sm_output = np.zeros_like(open_f.SM_SubX_m3_m3_value)
        sm_output = sm_output[:,:,0:7,:,:] #only want first 7 leads
        
        def find_percentile_of_score(sm_input, sm_output):
            for model in range(sm_input.shape[1]):
                for Y in range(sm_input.shape[3]):
                    for X in range(sm_input.shape[4]):
                        for idx,lead in enumerate(range(sm_input.shape[2])): 
                            #only looking at first 6 weeks for GMAO (only at weekly leads)
                            if idx == 0: #we can't take the weekly average of lead time 0
                                
                                all_values =all_rzsm_numpy[:,model,lead,Y,X]
                                # all_values = get_full_distribution(model,lead,Y,X)
                            
                                all_values = sorted(all_values)
                                all_values = [i for i in all_values if i !=0.0]
                            
                                #get percentile ranking (subtract 100 because higher anomalies are located in the upper percentiles)
                                sm_output[0,model,lead,Y,X] = pos(all_values,sm_input[0,model,lead,Y,X])
                                
                            elif idx%7 == 0:
                                #Find the weekly average
                                weekly_avg = np.mean(all_rzsm_numpy[:,model,idx-7:idx,Y,X],axis=1)
                                
                                all_values = sorted(weekly_avg)
                                all_values = [i for i in all_values if i !=0.0]
                                sm_output[0,model,idx//7,Y,X] = pos(all_values,sm_input[0,model,idx-7:idx,Y,X].mean())
                                
                        #if rzsm percentile is below the 20th: perform code else week is not in flash drought
            return(sm_output)
        
        sm_output = find_percentile_of_score(sm_input, sm_output)
        
        #Testing
        # model=3
        # lead=6
        # Y=17
        # X=44
        sm_input[0,3,:,17,44]
        sm_output[0,3,:,17,44]
        
        sm_input[0,3,:,10,10]
        sm_output[0,3,:,10,10]
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                RZSM_percentile = (['S','model','lead','Y','X'], sm_output[:,:,:,:,:]),
            ),
            coords = dict(
                X = open_f.X.values,
                Y = open_f.Y.values,
                lead = np.arange(0,7),
                S = open_f.S.values,
                model = open_f.model.values
            ),
            attrs = dict(
                Description = 'Soil Moisture RZSM percentiles from RZSM data.'),
        )                    
    

        var_final.to_netcdf(path = f'{percentile_dir}/rzsm_percentiles_{file[-14:-4]}.nc4')
        var_final.close()
        
        return()

    

#%%
'''Now that we have percentiles for RZSM anomalies, we can use the SMPD index
(Decrease from above the 40th percentile to below the 20th percentile)
'''

# os.chdir(smpd_dir)
# file_list = sorted(glob('anomaly_percentiles_*.nc4'))

# for file in file_list:
#     # file = file_list[0]
#     open_f = xr.open_dataset(file)

#     percentile_data = open_f.RZSM_percentile.to_numpy()
#     output_data = np.zeros_like(percentile_data)    
    
#     @njit(parallel=True)
#     def binary_occurence_smpd(percentile_data,output_data):
#         for model in prange(output_data.shape[1]):
#             for Y in range(output_data.shape[2]):
#                 for X in range(output_data.shape[4]):
#                     for lead in range(output_data.shape[2]):
#                     #if rzsm percentile is below the 20th: perform code else week is not in flash drought
#                         if percentile_data[0,model,lead,Y,X] <= 20:
#                         #if previous week is greater than 40th: week is in flash drought (1)
#                         #1 week code chunk (2 lines below)
#                             if percentile_data[0,model,lead-1,Y,X] >= 40:
#                                 output_data[0,model,lead,Y,X] = 1
#                         #must drop to below 20th percentile in 3 weeks or less
#                         #This is the 3-week code chunk below (next 5 lines)
#                             elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
#                                 percentile_data[0,model,lead-2,Y,X] <= 40 and \
#                                percentile_data[0,model,lead-3,Y,X] >= 40:
#                                     output_data[0,model,lead,Y,X]  = 1 #flash drought
#                         #2 weeks code chunk
#                             elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
#                                 percentile_data[0,model,lead-2,Y,X] >= 40:
#                                 output_data[0,model,lead,Y,X]   = 1 #flash drought   
#         return(output_data[:,:,:,:,:])
    
#     test = binary_occurence_smpd(percentile_data,output_data)

# np.count_nonzero(test == 1)

# percentile_data[0,0,:,11,37]



if __name__ == '__main__':
    p = Pool(num_processors)
    file_list_anomaly, all_rzsm_numpy, file_list_RZSM, all_rzsm_numpy = \
        get_files_and_data_setup(anom_dir=anom_dir,subx_RZSM_dir=subx_RZSM_dir)
    # p.map(create_anomaly_percentiles,file_list_anomaly)
    p.map(create_RZSM_percentiles,file_list_RZSM)