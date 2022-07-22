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
# num_processors = int(7)

dir1 = 'main_dir'
model_NAM1 = 'model_name'
num_processors = int('procs')

smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
os.chdir(smerge_dir)

SMERGE_RZSM_file = f'{smerge_dir}/smerge_sm_merged_remap.nc4'
SMERGE_anomaly_file = f'{smerge_dir}/RZSM_anomaly_SMERGE_merged.nc'

rzsm_file = xr.open_dataset(SMERGE_RZSM_file)
anomaly_file = xr.open_dataset(SMERGE_anomaly_file)

#%% GOOD
def create_anomaly_percentiles(anomaly_file,SMERGE_anomaly_file):
    
    try:
        outname = SMERGE_anomaly_file.split('/')[-1][:-3]
        xr.open_dataset(f'{outname}_percentiles.nc4')
        file_o = xr.open_dataset(f'{outname}_percentiles.nc4')
        return(file_o)
    except FileNotFoundError:
        print(f'Working on anomaly percentiles on file {SMERGE_anomaly_file}')      
        outname = SMERGE_anomaly_file.split('/')[-1]
        
        sm_anomaly = anomaly_file.CCI_ano.to_numpy()
        sm_output = np.zeros_like(sm_anomaly)
        
        '''Manually implement the source code for percentile of score (from scipy)
        Because numba does not like it. It is much faster with numba
        
        Source: https://github.com/scipy/scipy/blob/v1.8.1/scipy/stats/_stats_py.py#L1812-L1899'''
        
        @njit
        def find_percentile_of_score(sm_anomaly, sm_output):
            for Y in range(sm_anomaly.shape[1]):
                for X in range(sm_anomaly.shape[2]):
                    for day in range(sm_anomaly.shape[0]):
                        #Get all values for the specific X and Y for all days of the year
                        all_values =sm_anomaly[:,Y,X]
                        #remove np.nan
                        all_values = all_values[~np.isnan(all_values)]
                        
                        score = sm_anomaly[day,Y,X]
                        
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

                            sm_output[day,Y,X] = pct

            return(sm_output)
        
        sm_output = find_percentile_of_score(sm_anomaly, sm_output)

        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                Anomaly_percentile = (['time','Y','X'], sm_output[:,:,:]),
            ),
            coords = dict(
                X = anomaly_file.X.values,
                Y = anomaly_file.Y.values,
                time = anomaly_file.time.values,
            ),
            attrs = dict(
                Description = 'SMERGE Soil Moisture RZSM percentiles from anomaly data.'),
        )                    
    

        var_final.to_netcdf(path = f'{smerge_dir}/{outname}_percentiles.nc4')
        var_final.close()
        
        return()

#%% Good
def create_RZSM_m3_m3_percentiles(rzsm_file,SMERGE_RZSM_file):
    
    try:
        outname = SMERGE_RZSM_file.split('/')[-1][:-4]
        xr.open_dataset(f'{outname}_percentiles.nc4')
        file_o = xr.open_dataset(f'{outname}_percentiles.nc4')
        return(file_o)
    except FileNotFoundError:
        print(f'Working on anomaly percentiles on file {SMERGE_anomaly_file}')      
        #TODO: Perform a seven day-rolling mean first (has not yet been performed on data)
        rzsm = rzsm_file.rolling(time=7,center=True).mean()
        
        rzsm = rzsm.RZSM.to_numpy()
        rzsm_output = np.zeros_like(rzsm)
        
        '''Manually implement the source code for percentile of score (from scipy)
        Because numba does not like it. It is much faster with numba
        
        Source: https://github.com/scipy/scipy/blob/v1.8.1/scipy/stats/_stats_py.py#L1812-L1899'''
        
        @njit
        def find_percentile_of_score(rzsm, rzsm_output):
            for Y in range(rzsm_output.shape[1]):
                for X in range(rzsm_output.shape[2]):
                    for day in range(rzsm_output.shape[0]):
                        #Get all values for the specific X and Y for all days of the year
                        all_values =rzsm[:,Y,X]
                        #remove np.nan
                        all_values = all_values[~np.isnan(all_values)]
                        
                        score = rzsm[day,Y,X]
                        
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

                            rzsm_output[day,Y,X] = pct

            return(sm_output)
        
        sm_output = find_percentile_of_score(rzsm, rzsm_output)
        
     
        #Create new dataset
        var_final = xr.Dataset(
            data_vars = dict(
                RZSM_percentile = (['time','Y','X'], sm_output[:,:,:]),
            ),
            coords = dict(
                X = rzsm_file.X.values,
                Y = rzsm_file.Y.values,
                time = rzsm_file.time.values,
            ),
            attrs = dict(
                Description = 'SMERGE Soil Moisture RZSM percentiles from RZSM m3/m3 data.'),
        )                    
     
     
        var_final.to_netcdf(path = f'{smerge_dir}/{outname}_percentiles.nc4')
        var_final.close()
        
        return(var_final)


#%%


#%%
if __name__ == '__main__':
    anom_percentiles = create_anomaly_percentiles(anomaly_file,SMERGE_anomaly_file)
    rzsm_percentiles = create_RZSM_m3_m3_percentiles(rzsm_file,SMERGE_RZSM_file)
    '''We have the percentiles created from both RZSM m3/m3 and from the 42
    day anomalies. Next we need to classify flash drought with both of these files.
    '''
    
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
        print(f'Working on {file} for RZSM anomaly percentiles to SMPD binary classification.')

        @njit(parallel=True)
        def binary_occurence_smpd(percentile_data,output_data):
            for model in prange(output_data.shape[1]):
                for Y in range(output_data.shape[2]):
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
    # file = file_list[0]
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
        print(f'Working on {file} for RZSM m3/m3 percentiles to SMPD binary classification.')
        
        @njit(parallel=True)
        def binary_occurence_smpd(percentile_data,output_data):
            for model in prange(output_data.shape[1]):
                for Y in range(output_data.shape[2]):
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

if __name__ == '__main__':
    p = Pool(num_processors)
    file_list_anomaly, all_rzsm_numpy, file_list_RZSM, all_rzsm_numpy = \
        get_files_and_data_setup(anom_dir=anom_dir,subx_RZSM_dir=subx_RZSM_dir)
    # p.map(create_anomaly_percentiles,file_list_anomaly)
    p.map(create_RZSM_percentiles,file_list_RZSM)
    p.map(smpd_function_anomaly,file_list_anomaly_smpd)
    p.map(smpd_function_RZSM,file_list_rzsm_smpd)

