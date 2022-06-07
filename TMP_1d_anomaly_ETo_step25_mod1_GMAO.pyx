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
model_NUM = int('1')
model_NAM1 = 'GMAO'

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
def ETo_anomaly(int start_,int end_,int model_NUM,list init_date_list,str _date,str var):
    #_date=init_date_list[0]  
    
    # #Make sure that the models aren't overwriting each other 
    # if model_num == 1:
    #     os.system('sleep 2')
    # if model_num == 2:
    #     os.system('sleep 4')
    # if model_num == 3:
    #     os.system('sleep 6')
        
    #For each date, open each file and compute ETref with et
    #All files have the same initialized days (part of the pre-processing that is 
    #completed)
    def return_date_list():
        date_list = []
        for file in sorted(glob(f'{home_dir}/{var}*.nc4')):
            date_list.append(file[-14:-4])
        return(date_list)
            
    init_date_list = return_date_list()    

    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .npy in {home_dir}/{var}_anomaly_mod{model_NUM}.') 
    os.chdir(f'{home_dir}/{var}_anomaly_mod{model_NUM}')
    cdef int week_lead
    week_lead = 6

    #Used as the grid mask file to not iterate over useless grid cells
    eddi_file = xr.open_dataset(f'{home_dir}/EDDI_2011-06-14.nc4')
  
    #get julian day, timestamps, and the datetime
    def date_file_info(str _date,str variable):
        cdef int  start_julian, subtract, a_i
        cdef list a_date_out, a_julian_out,a_julian_out2
        
        open_f = xr.open_dataset(f'{home_dir}/{variable}_{_date}.nc4')
        a_date_in= open_f[f'{variable}'].lead.values
        #get the start date
        a_start_date = pd.to_datetime(open_f.S.values[0])
        #Add the dates based on index value in date_in
        a_date_out = []

        for a_i in range(len(a_date_in)):
            a_date_out.append(a_start_date + dt.timedelta(days = a_i))
        
        start_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
        #Julian day into a list                            
        a_julian_out = [start_julian + i for i in range(len(a_date_out))]
        
        if pd.to_datetime(_date).year % 4 == 0:
            subtract = 366
            a_julian_out2 = [i-subtract if i>366 else i for i in a_julian_out]
        else:
            subtract = 365
            a_julian_out2 = [i-subtract if i>365 else i for i in a_julian_out]
        
        #month of file
        INdate_for_month = dt.datetime(int(_date[0:4]),int(_date[5:7]),int(_date[8:10]))
        
        return(open_f,a_date_out,a_julian_out2,INdate_for_month,variable)
    
    subx,file_timestamp_list,file_julian_list,file_datetime,variable = date_file_info(_date=_date,variable='ETo')
     
    #Now convert to julian date and append coordinates
    subx2 = subx.assign_coords(lead = file_julian_list)
    
    # def moving_average(a, n=7) :
    #     ret = np.cumsum(a, dtype=float)
    #     ret[n:] = ret[n:] - ret[:-n]
    #     return ret[n - 1:] / n
    
    #Test on small subset
    #for i_Y in range(1):
    for i_Y in range(subx2[f'{variable}'].shape[3]):
        #for i_X in range(1):
         for i_X in range(subx2[f'{variable}'].shape[4]):
            if _date == '1999-01-10' and model_NUM==0 and start_ ==0:
                print(f'Working on lat {i_Y} and lon {i_X}')

            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE
            if (np.count_nonzero(np.isnan(eddi_file.EDDI[0,0,10,i_Y,i_X].values)) !=1)or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

                def dict1_subx2():
                    cdef dict summation_ETo_modN
                    cdef int idx,julian_d, end_julian, subtract, len_text, idx_lead
                    cdef list date_out, b_julian_out, b_julian_out2
                    
                    #append 7-day average from all files to a new dictionary, 
                    summation_ETo_modN = {}
                    
                    for idx,julian_d in enumerate(file_julian_list):
                        # if i==1:
                        #     break
                        #You must julian_d + week_lead because with ETo you need 7-day looking backwards into the past. Must have 7 values.
                        #Choose just model = 0 because we just need to know if there are 7-days total in any model
                        
                        
                        if idx % 7 == 0:
                            try:
                                if (len(subx2[f'{var}'].sel(lead=slice(julian_d,file_julian_list[idx+week_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)) == 7:
                                    summation_ETo_modN[f'{file_julian_list[idx+week_lead]}']=[]
                                    summation_ETo_modN[f'{file_julian_list[idx+week_lead]}'].append({f'{_date}':np.nanmean(subx2[f'{var}'].sel(lead=slice(julian_d,file_julian_list[idx+week_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)})   
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
                            
                            '''Convert lead dates to a vector, then add it back into a netcdf
                            because I cannot convert np.datetime to a pd.datetime, I will need
                            to iterate over each of the dates in a list from Et_ref'''
                            
                            #get the dates into a list
                            date_in= open_f[f'{variable}'].lead.values
                            #get the start date
                            start_date = pd.to_datetime(open_f.S.values[0])
                            #Add the dates based on index value in date_in
                            date_out = []
                            for i in range(len(date_in)):
                                date_out.append(start_date + dt.timedelta(days = i))
                            
                            #Convert to julian date
                            end_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
                            
                            b_julian_out = [end_julian + i for i in range(len(date_out))]
                            
                            '''Find out if that file has a leap year, subtract appropriately'''
                            if pd.to_datetime(file[-14:-4]).year % 4 == 0:
                                subtract = 366
                                b_julian_out2 = [i-subtract if i>366 else i for i in b_julian_out]
                            else:
                                subtract = 365
                                b_julian_out2 = [i-subtract if i>365 else i for i in b_julian_out]
                            
                        
                            Et_ref_open_f = open_f.assign_coords(lead = b_julian_out2)
                            
                            for idx,val in enumerate(b_julian_out2):
                                idx_lead = idx+week_lead
                               #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
                                if idx % 7 == 0:
                                    try:
                                        summation_ETo_modN[f'{b_julian_out2[idx_lead]}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.sel(lead=slice(val,b_julian_out2[idx_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)})   
                                   #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass
                    return(summation_ETo_modN)
                
                #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
                ETo_7_day_average_modN= dict1_subx2()
                                    
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
                    for idx,julian_date in enumerate(ETo_7_day_average_modN):
                        
                        out_eddi_dictionary[f'{julian_date}']=[]
                        
                        subset_by_date = ETo_7_day_average_modN[f'{julian_date}']
                        
                        mean_value = []
                        for date_init in range(len(subset_by_date)):
                            mean_value.append(list(subset_by_date[date_init].values())[0])
                            
                        mean_calc = np.nanmean(np.array(mean_value))
                        #calculate anomaly
                        anomaly_val = []
                        for date_init in range(len(subset_by_date)):
                            anomaly_val.append((list(subset_by_date[date_init].values())[0]-mean_calc))
                            
 
                        out_eddi_dictionary[f'{julian_date}'].append(anomaly_val)
                        
                    return(out_eddi_dictionary)
                    
                ETo_dict_modN = compute_anomaly(ETo_7_day_average_modN)

                '''Instead of re-looping through all of the files (very slow), we can start appending to 
                files one by one with the data that we have already collected.
                
                The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                has the initialized dates for each EDDI julian day summation'''
                
                def improve_RZSM_dictionary(dict ETo_7_day_average_modN, dict ETo_dict_modN):
                    cdef dict final_out_dictionary_all_eddi
                    cdef list sub_keys
                    cdef int idx, idxxx
                    cdef str julian_dattt
                    
                    final_out_dictionary_all_eddi = {}

                    for idx,julian_dattt in enumerate(ETo_7_day_average_modN):
                        final_out_dictionary_all_eddi[f'{julian_dattt}'] = [] #create an empty list to append to
                        sub_list = ETo_7_day_average_modN[f'{julian_dattt}'] #choose only the summation ETo with correct julian date
                        sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                        #sub list contains a dictionary for each julian date in current loop and the values of ETo
                        #Save the init date values for each julian date
                        for idxxx, init_date in enumerate(sub_list):
                            
                            sub_keys.append({list(init_date.keys())[0] :ETo_dict_modN[f'{julian_dattt}'][0][idxxx]})
                        
                        final_out_dictionary_all_eddi[f'{julian_dattt}'].append(sub_keys) 
                        
                    return(final_out_dictionary_all_eddi)
                
                ETo_next_dict_modN = improve_RZSM_dictionary(ETo_7_day_average_modN, ETo_dict_modN)

                
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_to_npy_file( ETo_next_dict_modN,  model_NUM):
                    cdef int idx_, index_val
                    cdef dict dic_init_and_eddi_val
                    cdef str init_day, i_val, fileOut
                    
                    for idx_,i_val in enumerate(ETo_next_dict_modN):
                      
                    #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict = ETo_next_dict_modN[f'{i_val}'][0]
                        
                        for dic_init_and_eddi_val in EDDI_final_dict:
                            # print(dic_init_and_eddi_val)
                            #Open up the file and insert the value
                            init_day = list(dic_init_and_eddi_val.keys())[0]
                            
                            '''When appending using several scripts, sometimes the file will
                            load before another file has finished filling in the file and causing
                            a ValueError: cannot reshape array of size 15328 into shape (2,4,45,27,59)'''
                            
                            lead_values = np.load(f'{var}_anomaly_{init_day}_julian_lead.npy',allow_pickle=True)
                            
                            #add values by julian day
                            # fileOut = f'{var}_anomaly_{init_day}.npy'
                            fileOut = ("{}_anomaly_{}.nc4".format(var,init_day))
                            
                            eddi_open = xr.open_dataset(fileOut)
                            index_val=np.where(lead_values == int(i_val))[0][0]
                            eddi_open.Variable[0,index_val,i_Y,i_X] = list(dic_init_and_eddi_val.values())[0]
                            
                            #For some reason it spits an error
   
                            eddi_open.to_netcdf(path = fileOut, mode ='a', engine='scipy')
                            eddi_open.close()
                            # del eddi_open
                            
                        
                            #OLD Way of doing it with .npy files
                            # fileOut = f'{var}_anomaly_{init_day}.npy'
                            
                            # try:
                            #     eddi_open = np.load(fileOut,allow_pickle=True)
                            # except ValueError:
                            #     os.system('sleep 2')
                            #     eddi_open = np.load(fileOut,allow_pickle=True)
                            
                            # lead_values = np.load(f'{var}_anomaly_{init_day}_julian_lead.npy',allow_pickle=True)
                            
                            # #Sometimes we found some indexes that really don't matter 
                            # try:
                            #     index_val=np.where(lead_values == int(i_val))[0][0]
                            #     eddi_open[0,index_val,i_Y,i_X] = list(dic_init_and_eddi_val.values())[0]
                            #     np.save(fileOut,eddi_open)
                            # except IndexError:
                            #     np.save(fileOut,eddi_open)
                            #     pass
                    
                add_to_npy_file(ETo_next_dict_modN,model_NUM)
                                
                #release some memory (doesn't work as well as I'd hoped)
                # del EDDI_next_dict_modN, EDDI_dict_modN, ETo_7_day_average_modN
                
    print(f'Completed date {_date} and saved into {home_dir}/{var}_anomaly_mod{model_NUM}.')
    #save the dates that were completed to not re-run
    if model_NUM == 3:
        os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_anomaly_npy_{model_NAM1}.txt')
    
    #Testing to see if this can help in releasing memory
    cdef str out = 'done'
    return(out)
#%%
'''For some odd reason, this function will keep allocating new memory (even though
there is nothing visible in the function that I can remove). To get around this,
it appears that I can only do 9 total dates between all 3 of my processes before
memory gets too high (potential memory leak), so add a break'''


#_date=init_date_list[0]
'''Read EDDI_completed_npy.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/{var}_completed_anomaly_npy_{model_NAM1}.txt',dtype='str')
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
        ETo_anomaly(start_, end_, model_NUM, init_date_list, _date,var)
        count+=1
        if count == 40:
            exit()


#Old way (this will help if I had to do every grid cell/lead anomalies)

# def run_loop(int count_total,int start_,int end_,int model_NUM,list init_date_list,str var):
#     '''Read ETo_completed_anomaly_npy.txt file to not have to re-run extra code'''
#     cdef int count
    
#     completed_dates = np.loadtxt(f'{script_dir}/ETo_completed_anomaly_npy_{model_NAM1}.txt',dtype='str')
#     try:
#         #first line contains a header, nothing with dates
#         completed_dates = completed_dates[:,1]
#     except IndexError:
#         completed_dates = ''
#     # completed_dates = pd.to_datetime(completed_dates[:],format='%Y-%m-%d')
#     #only work on dates that aren't completed
    
#     subset_completed_dates = [i[5:] for i in completed_dates]
    
#     #Save into a new directory after each date
#     #Make several new copies because it seems to break one of the files
#     new_directory = f'{home_dir}/{var}_anomaly_mod{model_NUM}/already_completed'

#     count=0
#     for _date in init_date_list[start_:end_]:    
#         if _date[5:] not in subset_completed_dates:
#             out_ = ETo_anomaly(start_, end_, model_NUM, init_date_list, _date,var)
#             #save
#             # os.system('sleep 30') #sometimes saving seems to not fully capture 1 file
#             # #Save multiple times because a file keeps getting broken
#             # for directory in [new_directory,new_directory2,new_directory3]:
#             #     os.system(f'cp {home_dir}/{var}_anomaly_mod{model_NUM}/*.nc4 {directory}/')
#             #     os.system('sleep 10')
                
#             break
#     return(count_total + 1)

# #Run in a seperate loop to possibly avoid memory issue, run_loop will add 1 to count_total
# #I tested count_total < 40 to keep the loop going, but it still eats up memory
# count_total = 0
# for i in range(2):
#     count_total = run_loop(count_total,start_,end_,model_NUM,init_date_list,var)
#     if i ==1:
#         #Break from script to see if memory leak is removed
#         gc.collect()
#         exit()

# # if __name__ == '__main__':
# #     call_function(start_, end_, model_NUM, init_date_list, _date)
