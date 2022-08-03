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
from numba import njit, prange


# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'
# var = 'RZSM'
# num_processors = int(2)

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'GMAO'
num_processors = int('7')

subx_RZSM_dir = f'{dir1}/Data/SubX/{model_NAM1}/SM_converted_m3_m3'

anom_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
percentile_dir = f'{anom_dir}/percentiles_anom_RZSM' 
smpd_dir = f'{anom_dir}/smpd_index'

os.system(f'mkdir -p {percentile_dir}')
os.system(f'mkdir -p {smpd_dir}')
var = 'RZSM'

#%% Create datasets for anomaly percentiles
def get_files_and_data_setup(anom_dir,subx_RZSM_dir):
    os.chdir(f'{anom_dir}') 
    
    global file_list_anomaly
    file_list_anomaly = sorted(glob('RZSM_anomaly*.nc4'))
    # file=file_list_anomaly[0]
    print('Loading RZSM anomaly files')
    all_anom = xr.open_mfdataset(file_list_anomaly, concat_dim=['S'], combine='nested')
    
    global all_anomaly_numpy
    all_anomaly_numpy = all_anom.RZSM_anom.to_numpy()
    
    global file_list_RZSM
    os.chdir(f'{subx_RZSM_dir}')
    file_list_RZSM = sorted(glob('SM_SubX_**.nc4'))
    
    print('Loading RZSM m3/m3 files')
    global all_rzsm
    all_rzsm = xr.open_mfdataset(file_list_RZSM, concat_dim=['S'], combine='nested')
    global all_rzsm_numpy 
    all_rzsm_numpy = all_rzsm.SM_SubX_m3_m3_value.to_numpy()
    
    return(file_list_anomaly,all_anomaly_numpy,file_list_RZSM,all_rzsm_numpy)
#%%

    
def create_anomaly_percentiles_SubX(file):
    
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
        
        
        @njit
        def find_percentile_of_score(sm_input, sm_output):
            for model in range(sm_input.shape[1]):
                for Y in range(sm_input.shape[3]):
                    for X in range(sm_input.shape[4]):
                        for lead in range(7): #only looking at first 6 weeks for GMAO
                            #Grab all the values within the distribution from all leads (about 13629 values)
                            all_values =all_anomaly_numpy[:,model,:,Y,X].flatten()
                            # all_values = get_full_distribution(model,lead,Y,X)
                            
                            all_values = sorted(all_values)
                            
                            #Remove 0.0. They shouldn't exist
                            all_values = [i for i in all_values if i !=0.0]
                                                    
                            score = sm_input[:,model,lead,Y,X]
                            
                            #construct percentiles with each day (percentile of score source code)
                            
                            a = np.asarray(all_values)
                            n = len(a)
                            if n == 0:
                                pass
                            else:
                                kind = 'rank' #from percentile of score function
                                if kind == 'rank':
                                    left = np.count_nonzero(a < score)
                                    right = np.count_nonzero(a <= score)
                                    pct = (right + left + (1 if right > left else 0)) * 50.0/n
    
                                    sm_output[0,model,lead,Y,X] = pct

            return(sm_output)

        sm_output_percentile = find_percentile_of_score(sm_input, sm_output)
        
    
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                Anomaly_percentile = (['S','model','lead','Y','X'], sm_output_percentile[:,:,0:7,:,:]),
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

    
#%% Create datasets for RZSM volumetric water (me/m3) content percentiles
def create_RZSM_percentiles_SubX(file):
    try:
        xr.open_dataset(f'{percentile_dir}/rzsm_percentiles_{file[-14:-4]}.nc4')
    except FileNotFoundError:
        print(f'Working on rzsm percentiles {file}')      
        # file=file_list_RZSM[0]
        #must change to directory to open file due to structure of list
        os.chdir(f'{subx_RZSM_dir}') 
        # file=file_list_anomaly[0]
        open_f = xr.open_dataset(file)
        
        sm_input = open_f.SM_SubX_m3_m3_value.to_numpy()
        sm_output = np.zeros_like(open_f.SM_SubX_m3_m3_value)
        sm_output = sm_output[:,:,0:7,:,:] #only want first 7 leads
        # all_rzsm_numpy = all_rzsm.SM_SubX_m3_m3_value.to_numpy()
        
        file_list_anomaly, all_anomaly_numpy, file_list_RZSM, all_rzsm_numpy ==get_files_and_data_setup(anom_dir,subx_RZSM_dir)

        
        '''Can't use numba for this unfortunately'''
        def find_percentile_of_score(sm_input, sm_output):
            for model in range(sm_input.shape[1]):
                for Y in range(sm_input.shape[3]):
                    for X in range(sm_input.shape[4]):
                        #only looking at first 6 weeks for GMAO (only at weekly leads)
                        
                        #Find the weekly average
                        weekly_avg1 = np.nanmean(all_rzsm_numpy[:,model,0:7,Y,X],axis=1)
                        weekly_avg2 = np.nanmean(all_rzsm_numpy[:,model,7:14,Y,X],axis=1)
                        weekly_avg3 = np.nanmean(all_rzsm_numpy[:,model,14:21,Y,X],axis=1)
                        weekly_avg4 = np.nanmean(all_rzsm_numpy[:,model,21:28,Y,X],axis=1)
                        weekly_avg5 = np.nanmean(all_rzsm_numpy[:,model,28:35,Y,X],axis=1)
                        weekly_avg6 = np.nanmean(all_rzsm_numpy[:,model,35:42,Y,X],axis=1)
                        
                        all_values_index_more_than_0 = np.concatenate([weekly_avg1,weekly_avg2,weekly_avg3,weekly_avg4,weekly_avg5, \
                                                     weekly_avg6])

                        all_values_index_more_than_0 = sorted(all_values_index_more_than_0)
                        all_values_index_more_than_0 = [i for i in all_values_index_more_than_0 if i !=0.0]
                        
                        for idx,lead in enumerate(range(sm_input.shape[2])): 

                            if idx == 0: #we can't take the weekly average of lead time 0
                                
                                all_values_index_0 =all_rzsm_numpy[:,model,0,Y,X] #get initial day values
                                
                                all_values_index_0 = np.concatenate([all_values_index_0,all_values_index_more_than_0])
                                # all_values = get_full_distribution(model,lead,Y,X)
                            
                                all_values_index_0 = sorted(all_values_index_0)
                                all_values_index_0 = [i for i in all_values_index_0 if i !=0.0]
                            
                                sm_output[0,model,idx,Y,X] = pos(all_values_index_0,sm_input[0,model,lead,Y,X])
                                
                            elif idx%7 == 0: #idx is the same as lead
                                #Find the weekly average

                                #Take weekly average as input into pos function
                                sm_output[0,model,idx//7,Y,X] = pos(all_values_index_more_than_0,sm_input[0,model,idx-7:idx,Y,X].mean())
                                
            return(sm_output)
        
        sm_output = find_percentile_of_score(sm_input, sm_output)
        
        #Testing
        # model=3
        # lead=6
        # Y=17
        # X=44
        # sm_input[0,3,:,17,44]
        # sm_output[0,3,:,17,44]
        
        # sm_input[0,3,:,10,10]
        # sm_output[0,3,:,10,10]
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

    

#%%Create SMPD binary occurence datasets for RZSM anomaly percentiles
'''Now that we have percentiles for RZSM anomalies, we can use the SMPD index
(Decrease from above the 40th percentile to below the 20th percentile)
'''

#test
# file = file_list_rzsm[0]
# for file in file_list_rzsm:
def smpd_function_anomaly(file):
    # file = file_list[0]
    os.chdir(percentile_dir)
    open_f = xr.open_dataset(file)

    percentile_data = open_f.anomaly_percentile.to_numpy()
    #file description (output)
    desc = 'SMPD index on RZSM anomaly data.'
    file_name = f'smpd_anomaly_{file[-14:-4]}.nc4'
    test_if_file_created = f'{smpd_dir}/{file_name}'
    
    output_data = np.zeros_like(percentile_data) 
    
    try:
        xr.open_dataset(test_if_file_created)
    except FileNotFoundError:
        print(f'Working on {file} to SMPD binary classification.')

        @njit(parallel=True)
        def binary_occurence_smpd(percentile_data,output_data):
            for model in prange(output_data.shape[1]):
                for Y in range(output_data.shape[3]):
                    for X in range(output_data.shape[4]):
                        for lead in range(output_data.shape[2]):
                        #if rzsm percentile is below the 20th: perform code else week is not in flash drought
                            if percentile_data[0,model,lead,Y,X] <= 20:
                            #if previous week is greater than 40th: week is in flash drought (1)
                            #1 week code chunk (2 lines below)
                                if percentile_data[0,model,lead-1,Y,X] >= 40:
                                    if lead != 0:
                                        #numpy will wrap around the end of the array,
                                        #this is messing up results
                                        output_data[0,model,lead,Y,X] = 1
                            #must drop to below 20th percentile in 3 weeks or less
                            #This is the 3-week code chunk below (next 5 lines)
                                elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-2,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-3,Y,X] >= 40:
                                        if lead != 1 and lead !=2 and lead !=3:
                                            output_data[0,model,lead,Y,X]  = 1 #flash drought
                            #2 weeks code chunk
                                elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-2,Y,X] >= 40:
                                        if lead != 2 and lead !=1:
                                            output_data[0,model,lead,Y,X]   = 1 #flash drought   
            return(output_data[:,:,:,:,:])
            
        # Y=0
        # X=14
        # percentile_data[0,model,:,Y,X] 
        # output_data[0,model,:,Y,X] 
        # OUTPUT[0,model,:,Y,X]  
        OUTPUT = binary_occurence_smpd(percentile_data,output_data)
    
        np.count_nonzero(OUTPUT == 1)
        
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                SMPD_RZSM = (['S','model','lead','Y','X'], OUTPUT[:,:,:,:,:]),
            ),
            coords = dict(
                X = open_f.X.values,
                Y = open_f.Y.values,
                lead = open_f.lead.values,
                S = open_f.S.values,
                model = open_f.model.values
            ),
            attrs = dict(
                Description = f'{desc}.'),
        )                    
    
    
        var_final.to_netcdf(path = f'{test_if_file_created}')
        var_final.close()
        
        return()
# percentile_data[0,0,:,11,37]

#  percentile_data[0,3,:,17,44]
#  output_data[0,3,:,17,44]
 
#  sm_input[0,3,:,10,10]
#  sm_output[0,3,:,10,10]

# Y=0
# X=13
# percentile_data[0,0,:,Y,X]
# percentile_data[0,0,:,0,14]
#%%Create SMPD binary occurence datasets for RZSM m3/m3 percentiles
def smpd_function_RZSM(file):
    # file = file_list_rzsm_smpd[0]
    open_f = xr.open_dataset(file)

    percentile_data = open_f.RZSM_percentile.to_numpy()
    output_data = np.zeros_like(percentile_data) 
    #file description (output)
    desc = 'SMPD index on RZSM m3/m3 data.'
    file_name = f'smpd_rzsm_{file[-14:-4]}.nc4'
    test_if_file_created = f'{smpd_dir}/{file_name}'
 
    try:
        xr.open_dataset(test_if_file_created)
    except FileNotFoundError:
        print(f'Working on {file} to SMPD binary classification.')
        
        @njit(parallel=True)
        def binary_occurence_smpd(percentile_data,output_data):
            for model in prange(output_data.shape[1]):
                for Y in range(output_data.shape[3]):
                    for X in range(output_data.shape[4]):
                        for lead in range(output_data.shape[2]):
                        #if rzsm percentile is below the 20th: perform code else week is not in flash drought
                            if percentile_data[0,model,lead,Y,X] <= 20:
                            #if previous week is greater than 40th: week is in flash drought (1)
                            #1 week code chunk (2 lines below)
                                if percentile_data[0,model,lead-1,Y,X] >= 40 and lead !=0:
                                    #numpy will wrap around the end of the array,
                                    #this is messing up results (lead !=0)
                                    output_data[0,model,lead,Y,X] = 1
                            #must drop to below 20th percentile in 3 weeks or less
                            #This is the 3-week code chunk below (next 5 lines)
                                elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-2,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-3,Y,X] >= 40:
                                        if lead != 1 and lead !=2 and lead !=3:
                                            output_data[0,model,lead,Y,X]  = 1 #flash drought
                            #2 weeks code chunk
                                elif percentile_data[0,model,lead-1,Y,X] <= 40 and \
                                    percentile_data[0,model,lead-2,Y,X] >= 40:
                                        if lead != 2 and lead !=1:
                                            output_data[0,model,lead,Y,X]   = 1 #flash drought   
            return(output_data[:,:,:,:,:])
            
        # Y=0
        # X=14
        # percentile_data[0,model,:,Y,X] 
        # output_data[0,model,:,Y,X] 
        # OUTPUT[0,model,:,Y,X]  
        OUTPUT = binary_occurence_smpd(percentile_data,output_data)
    
        np.count_nonzero(OUTPUT == 1)
        
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                SMPD_RZSM = (['S','model','lead','Y','X'], OUTPUT[:,:,:,:,:]),
            ),
            coords = dict(
                X = open_f.X.values,
                Y = open_f.Y.values,
                lead = open_f.lead.values,
                S = open_f.S.values,
                model = open_f.model.values
            ),
            attrs = dict(
                Description = f'{desc}.'),
        )                    
    
    
        var_final.to_netcdf(path = f'{test_if_file_created}')
        var_final.close()
            
        return('Complete')
#%%
#get files for functions
os.chdir(f'{percentile_dir}')
file_list_anomaly_smpd = sorted(glob('anomaly_percentiles_*.nc4'))
file_list_rzsm_smpd = sorted(glob('rzsm_percentiles_*.nc4'))


#Call all functions
if __name__ == '__main__':
    p = Pool(num_processors)
    global all_rzsm_numpy
    file_list_anomaly, all_anomaly_numpy, file_list_RZSM, all_rzsm_numpy = \
        get_files_and_data_setup(anom_dir=anom_dir,subx_RZSM_dir=subx_RZSM_dir)
    p.map(create_anomaly_percentiles_SubX,file_list_anomaly)
    p.map(create_RZSM_percentiles_SubX,file_list_RZSM)

    p.map(smpd_function_anomaly,file_list_anomaly_smpd)
    p.map(smpd_function_RZSM,file_list_rzsm_smpd)

#%%
# #Testing
# os.chdir('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/anomaly/smpd_index')

# file_list = os.listdir()
# total = 0
# file_occupied = []
# for file in file_list:
#     a=xr.open_dataset(file)
#     total+=a.SMPD_RZSM.sum()
#     if a.SMPD_RZSM.sum() > 200:
#         file_occupied.append(file)
        
# print(file_occupied)