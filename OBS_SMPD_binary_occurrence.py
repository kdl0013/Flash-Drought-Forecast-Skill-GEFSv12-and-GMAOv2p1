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
# num_processors = int(7)
dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
MERRA_dir = f'{dir1}/Data/MERRA2'

fileOUT_MERRA_root = f'{MERRA_dir}/SMPD_FD_MERRA.nc'
fileOUT_MERRA_root_anomaly = f'{MERRA_dir}/RZSM_anomaly_MERRA.nc'

CONUS_mask = xr.open_dataset(f'{dir1}/Scripts/CONUS_mask/NCA-LDAS_masks_SubX.nc4')
open_rzsm = xr.open_dataset(f'{MERRA_dir}/RZSM_anomaly_MERRA.nc')

anomaly_r = xr.zeros_like(open_rzsm).rename(RZSM_anom='SMPD') #OUTPUT
#Get datesopen_f
anomaly_date_r = pd.DataFrame(open_rzsm[list(open_rzsm.keys())[0]].time)
anomaly_date_r.index = pd.to_datetime(anomaly_date_r.iloc[:,0])
anomaly_date_r_list = list(anomaly_date_r.iloc[:,0])

def name(file):
    return(list(file.keys())[0])

#%%
def find_percentile_of_score(sm_input, sm_output):
    #For MERRA2 specifically 
    for Y in range(sm_input.shape[1]):
        print(Y)
        for X in range(sm_input.shape[2]):
            all_values =sm_input[:,Y,X]
            count=0
            for day in range(sm_input.shape[0]):

#%%

def find_percentile_of_score(sm_input, sm_output):
    #For MERRA2 specifically 
    for Y in range(sm_input.shape[1]):
        print(Y)
        for X in range(sm_input.shape[2]):
            all_values =sm_input[:,Y,X]
            count=0
            for day in range(sm_input.shape[0]):
                #(subtract from 100 because the anomalies are high)
                sm_output[day,Y,X] = pos(all_values,sm_input[day,Y,X])
                
                
                
                # print(pos(all_values,sm_input[day,Y,X]))
    return(sm_output)


sm_output = find_percentile_of_score(sm_input=open_rzsm[name(open_rzsm)].values, sm_output=anomaly_r[name(anomaly_r)].values)

#%%
def binary_occurence_smpd(array_):
    #Y-coords
    sm_output = np.zeros_like(array_)
    for i in range(array_.shape[1]):
        # X-coords
        for j in range(array_.shape[2]):
            for day in range(4,array_.shape[0]):
            #if rzsm percentile is below the 20th: perform code else week is not in flash drought
                if array_[day,i,j] <= 20:
                #if previous week is greater than 40th: week is in flash drought (1)
                #1 week code chunk (2 lines below)
                    if array_[day,i-1,j] >= 40:
                        sm_output[day,i,j] = 1
                #must drop to below 20th percentile in 3 weeks or less
                #This is the 3-week code chunk below (next 5 lines)
                    elif array_[day,i-1,j] <= 40 and array_[day,i-2,j] <= 40 and\
                        array_[day,i-3,j] >= 40:
                            sm_output[day,i,j] = 1 #flash drought
                #2 weeks code chunk
                    elif array_[day,i-1,j] <= 40 and array_[day,i-2,j] >= 40:
                        sm_output[day,i,j]  = 1 #flash drought   
    return(sm_output[:,:,:])
             
bin_output = binary_occurence_smpd(sm_output)

#%%
#make new arrays stacks arrays vertically
#want the binary_output (bin_output)s
#want a copy of the binary_output
#want a copy of the SM percentile values
def make_array(bin_output):
    calc_copy=bin_output.copy()
    arr_stack = np.vstack((calc_copy, sm_output ,calc_copy))  
    return(arr_stack)
  
caln=make_array(bin_output)

#Add weeks that are also in flash drought after the first week
#the algorithm ends when every grid cell is above the 20th percentile threshold
#repeat as necessary (until the output of each one is nearly the same)
#Add 1 week if the previous week is in flash drought and the current week is >= 0.54
@njit(parallel=True)
def add_weeks(array_):
    '''Split into seperate files because numba makes it difficult to pass in multiple items sometimes'''
    split_arr = int(len(array_)/3)
    binary_ = array_[:split_arr,]
    sm_val = array_[split_arr:-split_arr]
    binary_copy = array_[-split_arr:]
    for i in prange(binary_.shape[1]):
        for j in range(binary_.shape[2]):
            for day in range(binary_.shape[0]):
                if binary_[day,i,j]==0.0 and binary_[day,i-1,j]==1.0 and sm_val[day,i,j] <=20:
                        binary_copy[day,i,j] = 1  
    return(binary_copy[:,:,:])

#%%

#apply function
cal=add_weeks(caln)
print(sum(sum(cal[:,:,:])))

#complete a loop until there is no change in weeks that are also in flash drought (until diff =0)
count=0.0
iterations=0
diff=1
while diff!=0.0:
    caln=make_array(cal)
    #call function add_weeks
    cal=add_weeks(caln)
    print(sum(sum(cal[:,:])))
    cal_sum=sum(sum(cal[:,:]))
    diff=cal_sum-count
    count=cal_sum
    print(f'Difference between iterations is {diff}')
    iterations+=1
 
print(f'Number of iterations was {iterations}.')

@njit(parallel=True)
def remove_single_weeks(array_):
    cal_final=array_.copy()
    for j in prange(array_.shape[1]):
        for i in range(1,array_.shape[0]-1):
            if array_[i,j] == 1.0:
                if array_[i-1,j] == 0 and array_[i+1,j] == 0.0:
                    cal_final[i,j] = 0
    return(cal_final)

cal_out = remove_single_weeks(cal)
#%%

#Add weeks that are also in flash drought after the first week
#the algorithm ends when every grid cell is above the 20th percentile threshold
#repeat as necessary (until the output of each one is nearly the same)
#Add 1 week if the previous week is in flash drought and the current week is >= 0.54
@njit(parallel=True)
def add_weeks(array_):
    split_arr = int(len(array_)/3)
    binary_ = array_[:split_arr,]
    sm_val = array_[split_arr:-split_arr]
    binary_copy = array_[-split_arr:]
    for j in prange(binary_.shape[1]):
        for i in range(1,binary_.shape[0]):
            if binary_[i,j]==0.0 and binary_[i-1,j]==1.0 and sm_val[i,j] <=20:
                    binary_copy[i,j] = 1  
    return(binary_copy[:,:])

#%%
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
        sm_output = sm_output[:,:,0:7,:,:] #only want first 7 le
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

#%%
#Testing
os.chdir('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/anomaly/smpd_index')


 for file in file_list:
...     a=xr.open_dataset(file)
...     total+=a.SMPD_RZSM.sum()
...     if a.SMPD_RZSM.sum() > 0:
...             file_occupied.append(file)
