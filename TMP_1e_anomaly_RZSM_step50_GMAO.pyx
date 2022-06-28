#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly for RZSM.
    
Anomaly is calculated as the 7-day average of RZSM by 7-day window for each lead time.

E.g., First initialized file is 01-10-1999, and the actual first date for values is 01-11-1999
which has a julian day of 11. Find the weekly average of RZSM and find anomalies by week lead 
from all other years (total of 15 years (samples)).

Then subtract the mean from the value to find the anomaly.

This script only has to process 1 years worth of files because you need all files 
to create the distribution. While processing, it adds to all other files.

This script is broken up into different models with different step sizes for 
when the data should be calculated. Since running all files in parallel intially did work, 
it was very very slow because each file was being opened and re-opened many times and it
was creating a bottleneck.python

@author: kdl
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from scipy.stats import rankdata
import sys
import cython
import numpy
import gc



dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
start_ = int('50')
end_ = start_ + int('25')
model_NAM1 = 'GMAO'

# # Test for 1 step size and model
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# start_ = int('0')
# end_ = start_ + int('40')
# model_NUM = int('0')
# model_NAM1 = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
rzsm_dir = f'{home_dir}/SM_converted_m3_m3'

script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)
#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

gridMET_dir = f'{dir1}/Data/gridMET'

smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
###Files
file_list = os.listdir()

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    
global var
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


#%%    
'''process SubX files and create EDDI values'''
# def multiProcess_EDDI_SubX_TEST(_date):
def RZSM_anomaly(int start_, int end_,list init_date_list,str _date,str var):
    # _date=init_date_list[0]

    print(f'Calculating RZSM anomaly on SubX for {_date} and saving as .nc4 in {home_dir}.') 
    # os.chdir(f'{home_dir}/RZSM_anomaly_mod{model_NUM}')
    cdef int week_lead
    #Because of indexing, week lead actually needs to be 6 for slicing
    week_lead = 7

    #Used for eliminating iterating over grid cells that don't matter
    smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
    smerge_file_julian = smerge_file.copy()
    
    variable='SM_SubX_m3_m3'
    var_name = f'{variable}_value'
    
    subx = xr.open_dataset(f'{rzsm_dir}/{variable}_{_date}.nc4')
    
    # #get julian day, timestamps, and the datetime
    # def date_file_info( _date, variable):
    #     # cdef int start_julian, subtract
    #     # cdef list a_date_out, a_julian_out,a_julian_out2
        
    #     open_f = xr.open_dataset(f'{rzsm_dir}/{variable}_{_date}.nc4')
        
        
    #     a_date_in= open_f[f'{var_name}'].lead.values
    #     #get the start date
    #     a_start_date = pd.to_datetime(open_f.S.values[0])
    #     #Add the dates based on index value in date_in
    #     a_date_out = []

    #     for a_i in range(len(a_date_in)):
    #         a_date_out.append(a_start_date + dt.timedelta(days = a_i))
        
    #     start_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
    #     #Julian day into a list                            
    #     a_julian_out = [start_julian + i for i in range(len(a_date_out))]
        
    #     if pd.to_datetime(_date).year % 4 == 0:
    #         subtract = 366
    #         a_julian_out2 = [i-subtract if i>366 else i for i in a_julian_out]
    #     else:
    #         subtract = 365
    #         a_julian_out2 = [i-subtract if i>365 else i for i in a_julian_out]
        
    #     #month of file
    #     INdate_for_month = dt.datetime(int(_date[0:4]),int(_date[5:7]),int(_date[8:10]))
        
    #     return(open_f,a_date_out,a_julian_out2,INdate_for_month,var_name)
    
    # subx,file_timestamp_list,file_julian_list,file_datetime,var_name = date_file_info(_date=_date,variable='SM_SubX_m3_m3')
    
    # #Convert to julian day for processing
    # smerge_day = pd.to_datetime(smerge_file.CCI_ano.time.values)
    # smerge_julian = [i.timetuple().tm_yday for i in smerge_day]
    
    # smerge_file_julian= smerge_file_julian.assign_coords({'time': smerge_julian})
    
    
    #Now convert to julian date and append coordinates
    subx2 = subx.assign_coords(lead = np.arange(0,45))
    
    
    for i_Y in range(subx2[f'{var_name}'].shape[3]):
        for i_X in range(subx2[f'{var_name}'].shape[4]):
            if _date == '1999-01-10':
                print(f'Working on lat {i_Y} and lon {i_X}')

            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE
            
            if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1):
                
                def dict1_subx2():
                    # cdef dict summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3
                    # cdef int idx, end_julian, subtract, idx_lead
                    # cdef list date_out, b_julian_out, b_julian_out2, dates_to_keep
                    # cdef str file
                    
                    #append 7-day summation from all files to a new dictionary
                    #Technically, this is RZSM, but its easier to keep as ETo for EDDI, RZSM, and ETo
                    summation_ETo_mod0 = {}
                    summation_ETo_mod1 = {}
                    summation_ETo_mod2 = {}
                    summation_ETo_mod3 = {}
                    
                    for idx in range(subx2.lead.shape[0]):
                        #You must julian_d + week_lead because with RZSM you need 7-day looking backwards into the past. Must have 7 values.
                        #Choose just model = 0 because we just need to know if there are 7-days total in any model
                        # print(idx)
                        try:
                            if idx % 7 == 0 and idx !=0:
                                summation_ETo_mod0[f'{idx}']=[]
                                summation_ETo_mod0[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod1[f'{idx}']=[]
                                summation_ETo_mod1[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
    
                                summation_ETo_mod2[f'{idx}']=[]
                                summation_ETo_mod2[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
    
                                summation_ETo_mod3[f'{idx}']=[]
                                summation_ETo_mod3[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                        except IndexError:
                            pass
                        #Index error exists because there is no data at the end of the file
                        #That meets the current idx restrictions
        
                    '''7-day mean of RZSM by index:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''
                    
                    dates_to_keep = []
                    '''if within 3 months, keep files'''
                    for file in sorted(glob(f'{rzsm_dir}/SM*{_date[-5:]}.nc4')):
                            dates_to_keep.append(file)
                            
                    for file in dates_to_keep:
                        #Dont' re-open the same file
                        if file[-14:-4] != _date:
                            open_f = xr.open_dataset(file)
                            Et_ref_open_f = open_f.assign_coords(lead = np.arange(0,45))
                            
                            '''Now we need to append to the dictionary with the same julian date values'''
                            for idx in range(subx2.lead.shape[0]):
                                #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
                                if idx % 7 == 0 and idx != 0:
                                    try:
                                        summation_ETo_mod0[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod1[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod2[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod3[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.SM_SubX_m3_m3_value.isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   

                                    #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass
                            
                    return(summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3)
                
                #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
                SM_7_day_average_mod0,SM_7_day_average_mod1,SM_7_day_average_mod2,SM_7_day_average_mod3= dict1_subx2()
                
                
                                    
                '''Now we have created a dictionary that contains the:
                    1.) index value is the julian day
                    2.) list of dictionaries containing:
                    -init date file: summed 7-day value
                    [{'1999-01-10': 20.289343},  {'1999-01-15': 25.726818}, ....}]

                Now the dictionary is filled with all values from all files, we need to sort
                each julian day by having the largest values as number 1 ranking. Then we can append the 
                Et_out file with the proper julian day and the EDDI value
                
                '''
                
                '''Initially I thought I also needed to keep the mean for later reference, 
                but I do not'''
                def compute_anomaly(dict SM_7_day_average_modN):
                    cdef dict out_eddi_dictionary, out_mean_dictionary
                    cdef list mean_value, anomaly_val
                    cdef double mean_calc
                    cdef str julian_date
                    cdef int idx, date_init
                    
                    out_eddi_dictionary = {}
                    out_mean_dictionary = {}
                    #Key value in dictionary is julian_date
                    for idx,lead in enumerate(SM_7_day_average_modN):
                        
                        out_eddi_dictionary[f'{lead}']=[]
                        out_mean_dictionary[f'{lead}']=[]
                        
                        subset_by_date = SM_7_day_average_modN[f'{lead}']
                        
                        mean_value = []
                        for date_init in range(len(subset_by_date)):
                            mean_value.append(list(subset_by_date[date_init].values())[0])
                            
                        mean_calc = np.mean(np.array(mean_value))
                        #calculate anomaly
                        anomaly_val = []
                        for date_init in range(len(subset_by_date)):
                            anomaly_val.append((list(subset_by_date[date_init].values())[0]-mean_calc))
                            
 
                        out_eddi_dictionary[f'{lead}'].append(anomaly_val)
                        out_mean_dictionary[f'{lead}'].append(mean_calc)
                        
                    return(out_eddi_dictionary,out_mean_dictionary)
                    
                RZSM_dict_mod0, RZSM_mean_by_julian_day = compute_anomaly(SM_7_day_average_mod0)
                RZSM_dict_mod1, RZSM_mean_by_julian_day = compute_anomaly(SM_7_day_average_mod1)
                RZSM_dict_mod2, RZSM_mean_by_julian_day = compute_anomaly(SM_7_day_average_mod2)
                RZSM_dict_mod3, RZSM_mean_by_julian_day = compute_anomaly(SM_7_day_average_mod3)

                '''Instead of re-looping through all of the files (very slow), we can start appending to 
                files one by one with the data that we have already collected.
                
                The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                has the initialized dates for each EDDI julian day summation'''
                
                def improve_RZSM_dictionary(dict SM_7_day_average_modN,dict RZSM_dict_modN):
                    cdef dict final_out_dictionary_all_eddi
                    cdef list sub_keys
                    cdef int idx, idxxx
                    cdef str julian_dattt
                    
                    final_out_dictionary_all_eddi = {}

                    for idx,lead in enumerate(SM_7_day_average_modN):
                        final_out_dictionary_all_eddi[f'{lead}'] = [] #create an empty list to append to
                        sub_list = SM_7_day_average_modN[f'{lead}'] #choose only the summation ETo with correct julian date
                        sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                        #sub list contains a dictionary for each julian date in current loop and the values of ETo
                        #Save the init date values for each julian date
                        for idxxx, init_date in enumerate(sub_list):
                            
                            sub_keys.append({list(init_date.keys())[0] :RZSM_dict_modN[f'{lead}'][0][idxxx]})
                        
                        final_out_dictionary_all_eddi[f'{lead}'].append(sub_keys) 
                        
                    return(final_out_dictionary_all_eddi)
                
                RZSM_next_dict_mod0 = improve_RZSM_dictionary(SM_7_day_average_mod0, RZSM_dict_mod0)
                RZSM_next_dict_mod1 = improve_RZSM_dictionary(SM_7_day_average_mod1, RZSM_dict_mod1)
                RZSM_next_dict_mod2 = improve_RZSM_dictionary(SM_7_day_average_mod2, RZSM_dict_mod2)
                RZSM_next_dict_mod3 = improve_RZSM_dictionary(SM_7_day_average_mod3, RZSM_dict_mod3)

                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_to_nc_file(dict RZSM_next_dict_mod0,dict RZSM_next_dict_mod1,dict RZSM_next_dict_mod2,dict RZSM_next_dict_mod3):
                                        
                    for idx_,lead in enumerate(RZSM_next_dict_mod0):
                        lead = int(lead)
                    #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict0 = RZSM_next_dict_mod0[f'{lead}'][0]
                        EDDI_final_dict1 = RZSM_next_dict_mod1[f'{lead}'][0]
                        EDDI_final_dict2 = RZSM_next_dict_mod2[f'{lead}'][0]
                        EDDI_final_dict3 = RZSM_next_dict_mod3[f'{lead}'][0]
                        
                        for idx,dic_init_and_eddi_val in enumerate(EDDI_final_dict0):

                            # print(dic_init_and_eddi_val)
                            #Open up the file and insert the value
                            init_day = list(dic_init_and_eddi_val.keys())[0]
                            
                            var2 = 'RZSM'
                            #add values by julian day
                            fileOut = ("{}_anomaly_{}.nc".format(var2,init_day))
                            file_open = xr.open_dataset(fileOut)
                            file_open.close()

                            #Add data to netcdf file
                            file_open.RZSM_anom[0,0,lead,i_Y,i_X] = list(EDDI_final_dict0[idx].values())[0]
                            file_open.RZSM_anom[0,1,lead,i_Y,i_X] = list(EDDI_final_dict1[idx].values())[0]
                            file_open.RZSM_anom[0,2,lead,i_Y,i_X] = list(EDDI_final_dict2[idx].values())[0]
                            file_open.RZSM_anom[0,3,lead,i_Y,i_X] = list(EDDI_final_dict3[idx].values())[0]
   
                            file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                            file_open.close()
                            
                                           
                add_to_nc_file(RZSM_next_dict_mod0,RZSM_next_dict_mod1,RZSM_next_dict_mod2,RZSM_next_dict_mod3)
                
    print(f'Completed date {_date} and saved into {home_dir}.')
    
    #save the dates that were completed to not re-run
    os.system(f'echo Completed {_date} >> {script_dir}/RZSM_completed_anomaly_nc_{model_NAM1}.txt')
    
    return()
#%%



# _date=init_date_list[0]
'''Read RZSM_completed_anomaly_nc_.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/RZSM_completed_anomaly_nc_{model_NAM1}.txt',dtype='str')
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
        RZSM_anomaly(start_,end_,init_date_list, _date,var)



# count=0
# for _date in init_date_list[start_:end_]:    
#     if start_ == 50 and _date == '2000-01-10':
#         break
#     else:
#         RZSM_anomaly(start_,end_,init_date_list, _date,var)


