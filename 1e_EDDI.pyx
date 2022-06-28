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



dir1 = 'main_dir'
start_ = int('start_init')
end_ = start_ + int('init_step')
model_NAM1 = 'model_name'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# start_ = int('0')
# end_ = start_ + int('40')
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

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files
###Files
file_list = os.listdir()

global var
var = 'ETo'


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
def EDDI_function(int start_,int end_,list init_date_list,str _date,str var, HP_conus_mask):
    #_date=init_date_list[0]  
    variable='ETo' #Just used for creating dictionary from ETo files

    print(f'Calculating EDDI anomaly on SubX for {_date} and saving as .nc4 in {home_dir}.') 
  
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
    
    # subx,file_timestamp_list,file_julian_list,file_datetime,variable = date_file_info(_date=_date,variable='ETo')
     
    #Now convert to julian date and append coordinates
    subx = xr.open_dataset(f'{home_dir}/{variable}_{_date}.nc4')

    subx2 = subx.assign_coords(lead = np.arange(0,45))
    
    # def moving_average(a, n=7) :
    #     ret = np.cumsum(a, dtype=float)
    #     ret[n:] = ret[n:] - ret[:-n]
    #     return ret[n - 1:] / n
    
    #Test on small subset
    #for i_Y in range(1):
    for i_Y in range(subx2[f'{variable}'].shape[3]):
        #for i_X in range(1):
         for i_X in range(subx2[f'{variable}'].shape[4]):
            if _date == '1999-01-10':
                print(f'Working on lat {i_Y} and lon {i_X}')

            if HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7):
            #only work on grid cells with values like SMERGE
            # if (np.count_nonzero(np.isnan(subx2.EDDI[0,0,10,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

                def dict1_subx2():
                    cdef dict summation_ETo_mod0,summation_ETo_mod1,summation_ETo_mod2,summation_ETo_mod3
                    cdef int idx, end_julian, subtract, len_text, idx_lead
                    cdef list date_out, b_julian_out, b_julian_out2
                    
                    #append 7-day average from all files to a new dictionary
                    summation_ETo_mod0 = {}
                    summation_ETo_mod1 = {}
                    summation_ETo_mod2 = {}
                    summation_ETo_mod3 = {}
                    
                    for idx in range(subx2.lead.shape[0]):
                        #Only look at weekly leads
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
                               #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
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
                ETo_7_day_mod0,ETo_7_day_mod1,ETo_7_day_mod2,ETo_7_day_mod3= dict1_subx2()
                                    
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
                
                def compute_EDDI_on_ETo_7_day_sum_dict(dict ETo_7_day_modN):
                    cdef dict out_eddi_dictionary
                    cdef int idx
                    cdef float probability
                    cdef str lead
                    
                    out_eddi_dictionary = {}
                    #Key value in dictionary is julian_date
                    for idx,lead in enumerate(ETo_7_day_modN):
                        
                        out_eddi_dictionary[f'{lead}']=[]
                        
                        subset_by_date = ETo_7_day_modN[f'{lead}']
                        
                        #When looking at each julian date, now create a ranked list
                        list_rank = []
                        for date_init in range(len(subset_by_date)):
                            list_rank.append(list(subset_by_date[date_init].values())[0])
                        
                        ''' i = 1 for maximum ETo aggregation in time window'''
                        ranking = rankdata([-1 * i for i in list_rank]).astype(int)
                        tukey = (ranking - 0.33)/ (len(ranking) + 0.33)
                            
                        out_eddi = []
                        '''Now we can calculate EDDI'''
                        for probability in tukey:
                            if probability <= 0.5:
                                out_eddi.append(np.sqrt(-2*np.log(probability)))
                            else:
                                out_eddi.append((1-probability))  
                                    
                        #constants
                        c0 = 2.515517
                        c1 = 0.802853
                        c2 = 0.010328
                        d1 = 1.432788
                        d2 = 0.189269
                        d3 = 0.001308
                            
                        #Must reverse the sign of EDDI,so multiply by -1
                        final_out_eddi = []
                        for idx,w in enumerate(out_eddi):
                            final_out_eddi.append((w - ((c0 + c1*w + c2*w**2)/(1 + d1*w + d2*w**2 + d3*w**3)))*-1)
 
                        out_eddi_dictionary[f'{lead}'].append(final_out_eddi)
                        
                    return(out_eddi_dictionary)
                    
                EDDI_dict_mod0 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod0)
                EDDI_dict_mod1 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod1)
                EDDI_dict_mod2 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod2)
                EDDI_dict_mod3 = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_mod3)


                '''Instead of re-looping through all of the files (very slow), we can start appending to 
                files one by one with the data that we have already collected.
                
                The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                has the initialized dates for each EDDI julian day summation'''
                
                def improve_EDDI_dictionary(dict ETo_7_day_modN, dict EDDI_dict_modN):
                    cdef dict final_out_dictionary_all_eddi
                    cdef int idx,idxxx
                    cdef list sub_keys
                    
                    final_out_dictionary_all_eddi = {}

                    for idx,lead in enumerate(ETo_7_day_modN):
                        final_out_dictionary_all_eddi[f'{lead}'] = [] #create an empty list to append to
                        sub_list = ETo_7_day_modN[f'{lead}'] #choose only the summation ETo with correct julian date
                        sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                        #sub list contains a dictionary for each julian date in current loop and the values of ETo
                        #Save the init date values for each julian date
                        for idxxx, init_date in enumerate(sub_list):
                            
                            sub_keys.append({list(init_date.keys())[0] :EDDI_dict_modN[f'{lead}'][0][idxxx]})
                        
                        final_out_dictionary_all_eddi[f'{lead}'].append(sub_keys) 
                        
                    return(final_out_dictionary_all_eddi)
                
                EDDI_next_dict_mod0 = improve_EDDI_dictionary(ETo_7_day_mod0, EDDI_dict_mod0)
                EDDI_next_dict_mod1 = improve_EDDI_dictionary(ETo_7_day_mod1, EDDI_dict_mod1)
                EDDI_next_dict_mod2 = improve_EDDI_dictionary(ETo_7_day_mod2, EDDI_dict_mod2)
                EDDI_next_dict_mod3 = improve_EDDI_dictionary(ETo_7_day_mod3, EDDI_dict_mod3)

                
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_to_nc_file(dict EDDI_next_dict_mod0,dict EDDI_next_dict_mod1,dict EDDI_next_dict_mod2,dict EDDI_next_dict_mod3):
                    cdef int idx_, index_val
                    cdef dict dic_init_and_eddi_val
                    cdef str i_val
                    
                    for idx_,lead in enumerate(EDDI_next_dict_mod0):
                        lead = int(lead)

                        #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict0 = EDDI_next_dict_mod0[f'{lead}'][0]
                        EDDI_final_dict1 = EDDI_next_dict_mod1[f'{lead}'][0]
                        EDDI_final_dict2 = EDDI_next_dict_mod2[f'{lead}'][0]
                        EDDI_final_dict3 = EDDI_next_dict_mod3[f'{lead}'][0]
                        
                        for idx,dic_init_and_eddi_val in enumerate(EDDI_final_dict0):
                            # print(dic_init_and_eddi_val)
                            #Open up the file and insert the value
                            init_day = list(dic_init_and_eddi_val.keys())[0]
                            
                            fileOut = 'EDDI_{}.nc4'.format(init_day)
                            file_open = xr.open_dataset(fileOut)
                            file_open.close()

                            #Add data to netcdf file
                            file_open.EDDI[0,0,lead,i_Y,i_X] = list(EDDI_final_dict0[idx].values())[0]
                            file_open.EDDI[0,1,lead,i_Y,i_X] = list(EDDI_final_dict1[idx].values())[0]
                            file_open.EDDI[0,2,lead,i_Y,i_X] = list(EDDI_final_dict2[idx].values())[0]
                            file_open.EDDI[0,3,lead,i_Y,i_X] = list(EDDI_final_dict3[idx].values())[0]

                            file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                            file_open.close()
                          
                    
                add_to_nc_file(EDDI_next_dict_mod0,EDDI_next_dict_mod1,EDDI_next_dict_mod2,EDDI_next_dict_mod3)
                                
                #release some memory (doesn't work as well as I'd hoped)
                # del EDDI_next_dict_modN, EDDI_dict_modN, ETo_7_day_average_modN
                
    print(f'Completed date {_date} and saved into {home_dir}.')
    #save the dates that were completed to not re-run

    os.system(f'echo Completed {_date} >> {script_dir}/EDDI_completed_nc_{model_NAM1}.txt')

    return()
#%%
'''For some odd reason, this function will keep allocating new memory (even though
there is nothing visible in the function that I can remove). To get around this,
it appears that I can only do 9 total dates between all 3 of my processes before
memory gets too high (potential memory leak), so add a break'''


#_date=init_date_list[0]
'''Read EDDI_completed_npy.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/EDDI_completed_nc_{model_NAM1}.txt',dtype='str')

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
        EDDI_function(start_, end_, init_date_list, _date,var,HP_conus_mask)

