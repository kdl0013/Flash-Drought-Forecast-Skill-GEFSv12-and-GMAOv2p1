#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly mean for RZSM to later process for anomaly correlation
    coefficient (ACC).


@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob

dir1 = 'main_dir'
model_NAM1 = 'model_name'
start_ = int('start_init')
end_ = start_ + int('init_step')

# start_ = int('0')
# end_ = start_ + int('40')
# # Test for 1 step size and model
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
rzsm_dir = f'{home_dir}/SM_converted_m3_m3'

script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

new_dir_mean = 'mean_for_ACC'
os.system(f'mkdir {home_dir}/{new_dir_mean}')
#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

gridMET_dir = f'{dir1}/Data/gridMET'

smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
###Files
file_list = os.listdir()

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    
var='SM'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{rzsm_dir}/{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list()    

'''Steps for anomaly calculation. 
1.) For each date:
    2.) Go thorugh all files and find the 1 week average for the same dates:
        2a.) Find the mean
        3.) Calculate anomaly
        4.) Append all subsequent year's files since we have the first day 
    (only need to process the first year)
    

This code is re-factored from 1b_EDDI.py
    
'''

var='RZSM'
#%%    
'''process SubX files and create EDDI values'''
# def multiProcess_EDDI_SubX_TEST(_date):
def RZSM_mean(_date, var,anomaly_spread):
    # _date=init_date_list[0]
    print(f'Calculating RZSM mean for ACC on SubX for {_date} (and all subsequent years of the same date) and saving as .nc4 in {home_dir}/{new_dir_mean}.') 

    #Used for eliminating iterating over grid cells that don't matter
    smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
    smerge_file_julian = smerge_file.copy()
    
    variable='SM_SubX_m3_m3'
    var_name = f'{variable}_value'
    
    #Open all files for faster processing
    print('Currently we are not loading the dataset into memory because it eats up too much memory.')
    subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested')
    # subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested').persist()

    #Process indivdual files for output
    subx_out = xr.open_dataset(f'{home_dir}/{var}_SubX_m3_m3_{_date}.nc4')

    #Convert to julian day for processing
    smerge_day = pd.to_datetime(smerge_file.CCI_ano.time.values)
    smerge_julian = [i.timetuple().tm_yday for i in smerge_day]
    
    smerge_file_julian= smerge_file_julian.assign_coords({'time': smerge_julian})
     
    for i_Y in range(subx_out[f'{var_name}'].shape[3]):
        for i_X in range(subx_out[f'{var_name}'].shape[4]):
            if _date == '1999-01-10':
                print(f'Working on lat {i_Y} and lon {i_X}')

            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE
            
            if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1):
                
                def dict1_subx2(anomaly_spread):
                    
                    #append 7-day summation from all files to a new dictionary
                    #Technically, this is RZSM, but its easier to keep as ETo for EDDI, RZSM, and ETo
                    summation_ETo_mod0 = {}
                    summation_ETo_mod1 = {}
                    summation_ETo_mod2 = {}
                    summation_ETo_mod3 = {}
                    
                    for idx,julian_d in enumerate(subx_out.lead.values):
                        #You must julian_d + week_lead because with RZSM you need 7-day looking backwards into the past. Must have 7 values.
                        #Choose just model = 0 because we just need to know if there are 7-days total in any model
                        # print(idx)
                        try:
                            if idx % 7 == 0 and idx !=0:
                                summation_ETo_mod0[f'{julian_d}']=[]
                                summation_ETo_mod0[f'{julian_d}'].append({f'{_date}':np.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod1[f'{julian_d}']=[]
                                summation_ETo_mod1[f'{julian_d}'].append({f'{_date}':np.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
    
                                summation_ETo_mod2[f'{julian_d}']=[]
                                summation_ETo_mod2[f'{julian_d}'].append({f'{_date}':np.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
    
                                summation_ETo_mod3[f'{julian_d}']=[]
                                summation_ETo_mod3[f'{julian_d}'].append({f'{_date}':np.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                       
                        except IndexError:
                            pass
                        #Index error exists because there is no data at the end of the file
                        #That meets the current idx restrictions
        
                    '''7-day mean of RZSM by lead date for the first opened file:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''

                    dates_to_keep = []
                    '''grab yearly week lead files since we already have the distribution'''
                    for file in sorted(glob(f'{home_dir}/{var}_SubX*{_date[-5:]}.nc4')):
                            dates_to_keep.append(file)
                            
                    for file in dates_to_keep:
                        #Dont' re-open the same file
                        if file[-14:-4] != _date:
                            open_f = xr.open_dataset(file)

                            '''Now we need to append to the dictionary with the same julian date values'''
                            for idx,julian_d in enumerate(subx_out.lead.values):
                                #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
                                if idx % 7 == 0 and idx != 0:
                                    try:
                                        summation_ETo_mod0[f'{julian_d}'].append({f'{file[-14:-4]}':np.nanmean(open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod1[f'{julian_d}'].append({f'{file[-14:-4]}':np.nanmean(open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod2[f'{julian_d}'].append({f'{file[-14:-4]}':np.nanmean(open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod3[f'{julian_d}'].append({f'{file[-14:-4]}':np.nanmean(open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   

                                    #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass
                                    
                    def find_mean(summation_ETo_modN,mod_number,anomaly_spread):
                        out_dict = {}
                        '''This will return all the data for each file for averages'''
                        for k, julian_d in enumerate(summation_ETo_modN):

                            julian_d = int(julian_d)
                            if int(julian_d) <= anomaly_spread:
                                subtract_ = int(julian_d)-anomaly_spread
                                sub_small1 = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(366+subtract_,366))
                                sub_small2 = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(1,int(julian_d+anomaly_spread)))
                                sub_out = xr.concat([sub_small1,sub_small2],dim='lead')
                            elif julian_d >= (366-anomaly_spread):
                                add_ = 366-julian_d
                                diff_ = anomaly_spread - add_
                                sub_small1 = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(julian_d,julian_d+add_)) #get dates before new julian_day
                                sub_small2 = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(julian_d-anomaly_spread,julian_d))
                                sub_small3 = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(1,diff_))
                                sub_out = xr.concat([sub_small1,sub_small2,sub_small3],dim='lead')
                            else:
                                sub_out = subx_all[f'{var_name}'][:,mod_number,:,i_Y,i_X].sel(lead=slice(julian_d-anomaly_spread,julian_d+anomaly_spread)) 
 
                            mean_value = np.nanmean(sub_out)
                            out_dict[f'{julian_d}'] = mean_value
                            
                        
                            
                        return (out_dict)
                    
                    mean_mod0 = find_mean(summation_ETo_mod0,mod_number=0,anomaly_spread=anomaly_spread)
                    mean_mod1 = find_mean(summation_ETo_mod1,mod_number=1,anomaly_spread=anomaly_spread)
                    mean_mod2 = find_mean(summation_ETo_mod2,mod_number=2,anomaly_spread=anomaly_spread)
                    mean_mod3 = find_mean(summation_ETo_mod3,mod_number=3,anomaly_spread=anomaly_spread)
       
                    return(mean_mod0,mean_mod1,mean_mod2,mean_mod3,dates_to_keep)
                
                #Now we have the mean value for each year/week lead/model
                mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep= dict1_subx2(anomaly_spread=anomaly_spread)
                
               
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to RZSM'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found
                each mean_mod contains the mean values that need to be appended to all the years
                
                We are only keep 1 file since the mean is the same every year'''
                
                def add_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep):
                    
                    #get correct name of the file to be appended to 
                    
                    
                    for idx_,lead in enumerate(mean_mod0):
                        #We have the file, now open it
                        #add values by julian day
                        fileOut = "{}/{}_mean_{}".format(new_dir_mean,var,file_dates_to_keep[0][-14:])
                        file_open = xr.open_dataset(fileOut)
                        file_open.close()
                        
                    
                        #Add data to netcdf file
                        file_open.RZSM_mean[0,0,idx_*7,i_Y,i_X] = mean_mod0[f'{lead}']
                        file_open.RZSM_mean[0,1,idx_*7,i_Y,i_X] = mean_mod1[f'{lead}']
                        file_open.RZSM_mean[0,2,idx_*7,i_Y,i_X] = mean_mod2[f'{lead}']
                        file_open.RZSM_mean[0,3,idx_*7,i_Y,i_X] = mean_mod3[f'{lead}']

                add_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep)
                
    print(f'Completed date {_date} and saved into {home_dir}.')
    
    #save the dates that were completed to not re-run
    os.system(f'echo Completed {_date} >> {script_dir}/RZSM_completed_mean_nc_{model_NAM1}.txt')
    
    return()
#%%



# _date=init_date_list[0]
'''Read RZSM_completed_anomaly_nc_.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/RZSM_completed_mean_nc_{model_NAM1}.txt',dtype='str')
try:
    #first line contains a header, nothing with dates
    completed_dates = completed_dates[:,1]
except IndexError:
    completed_dates = ''
# completed_dates = pd.to_datetime(completed_dates[:],format='%Y-%m-%d')

#only work on dates that aren't completed
subset_completed_dates = [i[5:] for i in completed_dates]

for _date in init_date_list[start_:end_]:    
    if _date[5:] not in subset_completed_dates:
        RZSM_mean(_date=_date,var='RZSM',anomaly_spread=42)



# count=0
# for _date in init_date_list[start_:end_]:    
#     if start_ == 50 and _date == '2000-01-10':
#         break
#     else:
#         RZSM_anomaly(start_,end_,init_date_list, _date,var)


