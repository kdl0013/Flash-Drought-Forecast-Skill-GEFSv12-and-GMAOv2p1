#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly for RZSM and EDDI

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
start_ = int('25')
end_ = start_ + int('25')
model_NAM1 = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# start_ = int('0')
# end_ = start_ + int('25')
# model_NUM = int('0')
# model_NAM1 = 'GMAO'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
script_dir = f'{dir1}/Scripts'

os.chdir(home_dir)
#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'
gridMET_dir = f'{dir1}/Data/gridMET'
smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
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
def ETo_anomaly(int start_,int end_,list init_date_list,str _date,str var):
    #_date=init_date_list[0]  
  
    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .nc4 in {home_dir}') 
    #Because of indexing, week lead actually needs to be 6 for slicing
    cdef int week_lead
    week_lead = 7

    #Used as the grid mask file to not iterate over useless grid cells
    #Used for eliminating iterating over grid cells that don't matter
    smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
    smerge_file_julian = smerge_file.copy()
    
    variable='ETo'
    subx = xr.open_dataset(f'{home_dir}/{variable}_{_date}.nc4')
  
    # #get julian day, timestamps, and the datetime
    # def date_file_info(str _date,str variable):
    #     cdef int  start_julian, subtract, a_i
    #     cdef list a_date_out, a_julian_out,a_julian_out2
        
    #     open_f = xr.open_dataset(f'{home_dir}/{variable}_{_date}.nc4')
    #     a_date_in= open_f[f'{variable}'].lead.values
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
        
    #     return(open_f,a_date_out,a_julian_out2,INdate_for_month,variable)
    
    # subx,file_timestamp_list,file_julian_list,file_datetime,var_name = date_file_info(_date=_date,variable='ETo')
     
    #Now convert to julian date and append coordinates
    subx2 = subx.assign_coords(lead = np.arange(0,45))
    
    #Test on small subset
    #for i_Y in range(1):
    for i_Y in range(subx2[f'{variable}'].shape[3]):
        #for i_X in range(1):
         for i_X in range(subx2[f'{variable}'].shape[4]):
            if _date == '1999-01-10':
                print(f'Working on lat {i_Y} and lon {i_X}')

            #only work on grid cells with values like SMERGE
            if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

                def dict1_subx2():
                    cdef dict summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3
                    cdef int idx, len_text
                    
                    #append 7-day average from all files to a new dictionary
                    summation_ETo_mod0 = {}
                    summation_ETo_mod1 = {}
                    summation_ETo_mod2 = {}
                    summation_ETo_mod3 = {}

                    
                    for idx in range(subx2.lead.shape[0]):

                        if idx % 7 == 0 and idx !=0:
                            try:
                                summation_ETo_mod0[f'{idx}']=[]
                                summation_ETo_mod0[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{variable}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod1[f'{idx}']=[]
                                summation_ETo_mod1[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{variable}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod2[f'{idx}']=[]
                                summation_ETo_mod2[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{variable}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod3[f'{idx}']=[]
                                summation_ETo_mod3[f'{idx}'].append({f'{_date}':np.nanmean(subx2[f'{variable}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})
                            except IndexError:
                                pass

                    '''7-day mean of EDDI by index:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''
                    
                    dates_to_keep = []
                    '''if within 3 months, keep files'''
                    len_text=18
                    #ETo_anomaly is also in this directory, we don't want those
                    for file in sorted(glob(f'{home_dir}/ETo*{_date[-5:]}.nc4')):
                        if len(file.split('/')[-1]) == len_text:
                            dates_to_keep.append(file)
                                        
                    for file in dates_to_keep:
                        #Dont' re-open the same file
                        if file[-14:-4] != _date:
                            #Open up ETo file
                            open_f = xr.open_dataset(file)                        
                            Et_ref_open_f = open_f.assign_coords(lead = np.arange(0,45))
                            
                            for idx in range(subx2.lead.shape[0]):
                                
                                if idx % 7 == 0 and idx != 0:
                                    try:
                                        summation_ETo_mod0[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod1[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod2[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod3[f'{idx}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                                        
                                   #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass
                    return(summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3)
                
                #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
                
                ETo_7_day_average_mod0,ETo_7_day_average_mod1,ETo_7_day_average_mod2,ETo_7_day_average_mod3= dict1_subx2()
                                    
                '''Now we have created a dictionary that contains the:
                    1.) index value is the julian day
                    2.) list of dictionaries containing:
                    -init date file: summed 7-day value
                    [{'1999-01-10': 20.289343},  {'1999-01-15': 25.726818}, ....}]

                Now the dictionary is filled with all values from all files, we need to sort
                each julian day by having the largest values as number 1 ranking. Then we can append the 
                Et_out file with the proper julian day and the EDDI value
                
                EDDI calculation source :
                    M. Hobbins, A. Wood, D. McEvoy, J. Huntington, C. Morton, M. Anderson, and C. Hain (June 2016): The Evaporative Demand Drought Index: Part I – Linking Drought Evolution to Variations in Evaporative Demand. J. Hydrometeor., 17(6),1745-1761, doi:10.1175/JHM-D-15-0121.1.
                '''
                def compute_anomaly(dict ETo_7_day_average_modN):
                    cdef dict out_eddi_dictionary
                    cdef list mean_value, anomaly_val
                    cdef float mean_calc
                    cdef int idx, date_init
                    cdef str julian_date
                    
                    out_eddi_dictionary = {}
                    #Key value in dictionary is julian_date
                    for idx,lead in enumerate(ETo_7_day_average_modN):
                        
                        out_eddi_dictionary[f'{lead}']=[]
                        
                        subset_by_date = ETo_7_day_average_modN[f'{lead}']
                        
                        mean_value = []
                        for date_init in range(len(subset_by_date)):
                            mean_value.append(list(subset_by_date[date_init].values())[0])
                            
                        mean_calc = np.nanmean(np.array(mean_value))
                        #calculate anomaly
                        anomaly_val = []
                        for date_init in range(len(subset_by_date)):
                            anomaly_val.append((list(subset_by_date[date_init].values())[0]-mean_calc))
                            
 
                        out_eddi_dictionary[f'{lead}'].append(anomaly_val)
                        
                    return(out_eddi_dictionary)
                    
                ETo_dict_mod0 = compute_anomaly(ETo_7_day_average_mod0)
                ETo_dict_mod1 = compute_anomaly(ETo_7_day_average_mod1)
                ETo_dict_mod2 = compute_anomaly(ETo_7_day_average_mod2)
                ETo_dict_mod3 = compute_anomaly(ETo_7_day_average_mod3)


                '''Instead of re-looping through all of the files (very slow), we can start appending to 
                files one by one with the data that we have already collected.
                
                The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                has the initialized dates for each EDDI julian day summation'''
                
                def improve_RZSM_dictionary(dict ETo_7_day_average_modN, dict ETo_dict_modN):
                    cdef dict final_out_dictionary_all_eddi
                    cdef list sub_keys
                    cdef int idx, idxxx

                    
                    final_out_dictionary_all_eddi = {}

                    for idx,lead in enumerate(ETo_7_day_average_modN):
                        final_out_dictionary_all_eddi[f'{lead}'] = [] #create an empty list to append to
                        sub_list = ETo_7_day_average_modN[f'{lead}'] #choose only the summation ETo with correct julian date
                        sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                        #sub list contains a dictionary for each julian date in current loop and the values of ETo
                        #Save the init date values for each julian date
                        for idxxx, init_date in enumerate(sub_list):
                            
                            sub_keys.append({list(init_date.keys())[0] :ETo_dict_modN[f'{lead}'][0][idxxx]})
                        
                        final_out_dictionary_all_eddi[f'{lead}'].append(sub_keys) 
                        
                    return(final_out_dictionary_all_eddi)
                
                ETo_next_dict_mod0 = improve_RZSM_dictionary(ETo_7_day_average_mod0, ETo_dict_mod0)
                ETo_next_dict_mod1 = improve_RZSM_dictionary(ETo_7_day_average_mod1, ETo_dict_mod1)
                ETo_next_dict_mod2 = improve_RZSM_dictionary(ETo_7_day_average_mod2, ETo_dict_mod2)
                ETo_next_dict_mod3 = improve_RZSM_dictionary(ETo_7_day_average_mod3, ETo_dict_mod3)


                
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_to_nc_file(dict ETo_next_dict_mod0,dict ETo_next_dict_mod1,dict ETo_next_dict_mod2,dict ETo_next_dict_mod3):
                    cdef int idx_, index_val
                    cdef dict dic_init_and_eddi_val
                    cdef str init_day, i_val, fileOut
                    
                    for idx_,lead in enumerate(ETo_next_dict_mod0):
                        lead = int(lead)
                      
                    #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict0 = ETo_next_dict_mod0[f'{lead}'][0]
                        EDDI_final_dict1 = ETo_next_dict_mod1[f'{lead}'][0]
                        EDDI_final_dict2 = ETo_next_dict_mod2[f'{lead}'][0]
                        EDDI_final_dict3 = ETo_next_dict_mod3[f'{lead}'][0]

                        
                        for idx,dic_init_and_eddi_val in enumerate(EDDI_final_dict0):
                            # print(dic_init_and_eddi_val)
                            #Open up the file and insert the value
                            init_day = list(dic_init_and_eddi_val.keys())[0]

                            #add values by julian day
                            fileOut = ("{}_anomaly_{}.nc".format(variable,init_day))
                            file_open = xr.open_dataset(fileOut)
                            file_open.close()
    
                            #Add data to netcdf file
                            file_open.ETo_anom[0,0,lead,i_Y,i_X] = list(EDDI_final_dict0[idx].values())[0]
                            file_open.ETo_anom[0,1,lead,i_Y,i_X] = list(EDDI_final_dict1[idx].values())[0]
                            file_open.ETo_anom[0,2,lead,i_Y,i_X] = list(EDDI_final_dict2[idx].values())[0]
                            file_open.ETo_anom[0,3,lead,i_Y,i_X] = list(EDDI_final_dict3[idx].values())[0]
   
                            file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                            file_open.close()
                          
                    
                add_to_nc_file(ETo_next_dict_mod0,ETo_next_dict_mod1,ETo_next_dict_mod2,ETo_next_dict_mod3)
                
                
    #save the dates that were completed to not re-run
    os.system(f'echo Completed {_date} >> {script_dir}/ETo_completed_anomaly_nc_{model_NAM1}.txt')
    print(f'Completed date {_date} and saved into {home_dir}.')
           
    # os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_anomaly_nc_{model_NAM1}.txt')
    return()
#END FUNCTION
#%%
#Call function

#_date=init_date_list[0]
'''Read ETo_completed_npy.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/{var}_completed_anomaly_nc_{model_NAM1}.txt',dtype='str')

try:
    #first line contains a header, nothing with dates
    completed_dates = completed_dates[:,1]
except IndexError:
    completed_dates = ''
# completed_dates = pd.to_datetime(completed_dates[:],format='%Y-%m-%d')
#only work on dates that aren't completed

subset_completed_dates = [i[5:] for i in completed_dates]

count=0
for _date in init_date_list[start_:end_]:
    if _date[5:] not in subset_completed_dates:
        ETo_anomaly(start_, end_, init_date_list, _date,var)

    