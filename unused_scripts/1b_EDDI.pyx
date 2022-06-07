#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate EDDI using Hobbins et al. (2016) method
https://journals.ametsoc.org/view/journals/hydr/17/6/jhm-d-15-0121_1.xml    

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
from setuptools import setup, Extension
import numpy
from Cython.Distutils import build_ext as _build_ext


dir1 = 'main_dir'
start_ = int('start_init')
end_ = start_ + int('init_step')
model_NUM = int('model_number')
model_NAM1 = 'model_name'

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

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)

'''Steps for EDDI calculation. 
1.) For each date:
    2.) Go thorugh all files and find the 1 week summed ETo for the same dates:
        2a.) Only append the first file with rank and EDDI
        3.) Then rank them according to which weeks have the highest summation:
        4.) Calculate Tukey plotting position
        5.) Calculate EDDI
        6.) Append all subsequent year's files since we have the first day 
    (only need to process the first year)
    
    
To make the multiprocessing function faster (because it appears that accessing the same file is
actually slowing down the processing), I am creating many instances of the 
ETo_ SubX files into new directories. Then I can np.random a number and then
we can pull from that directory.
    
'''
global var
var = 'ETo'

#Same function inside function below, but need to keep a copy of it outside
def return_date_list():
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list()   
#%%    
'''process SubX files and create EDDI values'''
# def multiProcess_EDDI_SubX_TEST(_date):
def EDDI_function(str start_,str end_,int model_NUM,list init_date_list,str _date,str var):
    #_date=init_date_list[0]
     
    def return_date_list():
        date_list = []
        for file in sorted(glob(f'{var}*.nc4')):
            date_list.append(file[-14:-4])
        return(date_list)
            
    init_date_list = return_date_list()    
    

    print(f'Calculating EDDI on SubX for {_date} and saving as .npy in {home_dir}/EDDI_mod{model_NUM}.') 
    os.chdir(f'{home_dir}/EDDI_mod{model_NUM}')
    
    cdef int week_lead
    week_lead = 6

    #Used for eliminating iterating over grid cells that don't matter
    smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
  

    #get julian day, timestamps, and the datetime
    def date_file_info(str _date,str variable):
        cdef int start_julian, subtract
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
        
        return(open_f,a_date_out,a_julian_out2,INdate_for_month)
    
    subx,file_timestamp_list,file_julian_list,file_datetime = date_file_info(_date=_date,variable='ETo')
     
    #Now convert to julian date and append coordinates
    subx2 = subx.assign_coords(lead = file_julian_list)
    
    #Test on small subset
    #for i_Y in range(1):
    for i_Y in range(subx2.ETo.shape[3]):
        #for i_X in range(1):
         for i_X in range(subx2.ETo.shape[4]):
            if _date == '1999-01-10' and model_NUM==0 and start_ ==0:
                print(f'Working on lat {i_Y} and lon {i_X}')

                          
            '''Next if statement creates a mask to restrict within CONUS, because
            EDDI covers North America, this will create the mask using SMERGE plus 2 additional grid cell.
            Don't loop over the values that shouldn't have values.
            SMERGE indexes 38,6 and 38,7 are missing, but we can compute EDDI values 
            
            if below statement is FALSE, the file value is already np.nan '''
            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            if (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7)):

                def dict1_subx2():
                    cdef dict summation_ETo_modN
                    cdef int idx,julian_d, end_julian, subtract
                    cdef list date_out, b_julian_out, b_julian_out2
                    
                    #append 7-day summation from all files to a new dictionary, 
                    summation_ETo_modN = {}
                    
                    for idx,julian_d in enumerate(file_julian_list):
                        #You must julian_d + week_lead because with EDDI you need 7-day looking backwards into the past. Must have 7 values.
                        #Choose just model = 0 because we just need to know if there are 7-days total in any model
                        
                        if idx < 39:
                            if (len(subx2.ETo.sel(lead=slice(julian_d,file_julian_list[idx+week_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)) == 7:
                                summation_ETo_modN[f'{file_julian_list[idx+week_lead]}']=[]
                                summation_ETo_modN[f'{file_julian_list[idx+week_lead]}'].append({f'{_date}':np.nansum(subx2.ETo.sel(lead=slice(julian_d,file_julian_list[idx+week_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)})   

                    '''7-day sum of Eto by index:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''
                    
                    #Keep only files within a certain date
                    def limit_file_list(str _date):
                        cdef int month_start, test_1, test_2
                        cdef list dates_to_keep
                        cdef str file
                        
                        month_start = pd.to_datetime(_date).month
                        
                        dates_to_keep = []
                        '''if within 3 months, keep files'''
                        for file in sorted(glob(f'{home_dir}/ETo*.nc4')):
                            test_1 = month_start - pd.to_datetime(file[-14:-4]).month
                            test_2 = pd.to_datetime(file[-14:-4]).month - month_start
                            
                            if test_1 <=-9 or test_1>= 9:
                                dates_to_keep.append(file)
                            elif abs(test_1) <=3 or abs(test_2) <=3:
                                dates_to_keep.append(file)
                                
                        return (dates_to_keep)
                    
                    '''DON'T USE THIS INCOMPLETE FUNCTION, CPU USAGE DOES NOT ALLOW FOR ME TO RUN EXTRA FILES
                    (WHICH THIS FUNCTION COULD ALLOW ME TO DO)
                    This function could work but doesn't really account for the overlapping julian days between years
                    #Keep only files within a certain date
                    def limit_file_list_days(_date):
                        day_of_year = pd.to_datetime(_date).timetuple().tm_yday
                        
                        #Number of days between files
                        #change file year to same year
                        new_file = _date[0:4] + file[4:]
                        num_days = abs(pd.to_datetime(_date) - pd.to_datetime(new_file)).days
                        
                         num_days %365
                        
                        dates_to_keep = []
                        #if within 60 days, keep files
                        file = sorted(glob(f'{home_dir}/ETo*.nc4'))[20]
                        for file in sorted(glob(f'{home_dir}/ETo*.nc4')):
                            test_1 = day_of_year - pd.to_datetime(file[-14:-4]).timetuple().tm_yday
                            
                            if abs(test_1) <=60:
                                dates_to_keep.append(file)
                                
                        return (dates_to_keep)
                    '''
                    
                    dates_to_keep = limit_file_list(_date)
                    
                    for file in dates_to_keep:
                        #Dont' re-open the same file
                        if file[-14:-4] != _date:
                            
                            open_f = xr.open_dataset(file)
                            
                            '''Convert lead dates to a vector, then add it back into a netcdf
                            because I cannot convert np.datetime to a pd.datetime, I will need
                            to iterate over each of the dates in a list from Et_ref'''
                            
                            #get the dates into a list
                            date_in= open_f.ETo.lead.values
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
                            
                            '''Now we need to append to the dictionary with the same julian date values'''
                            for idx,val in enumerate(b_julian_out2):
                                #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
                                if idx < 39:
                                    try:
                                        summation_ETo_modN[f'{b_julian_out2[idx+week_lead]}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,b_julian_out2[idx+week_lead])).isel(S=0, model=model_NUM, X=i_X, Y=i_Y).values)})   
                                    except KeyError:
                                        pass
                            
                    return(summation_ETo_modN)
                
                #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
                ETo_7_day_modN= dict1_subx2()
                                    
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
                    cdef list mean_value, anomaly_val, final_out_eddi
                    cdef double ranking, tukey, c0, c1, c2, d1, d2, d3
                    
                    out_eddi_dictionary = {}
                    #Key value in dictionary is julian_date
                    for idx,julian_date in enumerate(ETo_7_day_modN):
                        
                        out_eddi_dictionary[f'{julian_date}']=[]
                        
                        subset_by_date = ETo_7_day_modN[f'{julian_date}']
                        
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
 
                        out_eddi_dictionary[f'{julian_date}'].append(final_out_eddi)
                        
                    return(out_eddi_dictionary)
                    
                EDDI_dict_modN = compute_EDDI_on_ETo_7_day_sum_dict(ETo_7_day_modN)

                '''Instead of re-looping through all of the files (very slow), we can start appending to 
                files one by one with the data that we have already collected.
                
                The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                has the initialized dates for each EDDI julian day summation'''
                
                def improve_EDDI_dictionary(dict ETo_7_day_modN, dict EDDI_dict_modN):
                    cdef dict final_out_dictionary_all_eddi
                    cdef list sub_keys
                    
                    final_out_dictionary_all_eddi = {}

                    for idx,julian_dattt in enumerate(ETo_7_day_modN):
                        final_out_dictionary_all_eddi[f'{julian_dattt}'] = [] #create an empty list to append to
                        sub_list = ETo_7_day_modN[f'{julian_dattt}'] #choose only the summation ETo with correct julian date
                        sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                        #sub list contains a dictionary for each julian date in current loop and the values of ETo
                        #Save the init date values for each julian date
                        for idxxx, init_date in enumerate(sub_list):
                            
                            sub_keys.append({list(init_date.keys())[0] :EDDI_dict_modN[f'{julian_dattt}'][0][idxxx]})
                        
                        final_out_dictionary_all_eddi[f'{julian_dattt}'].append(sub_keys) 
                        
                    return(final_out_dictionary_all_eddi)
                
                EDDI_next_dict_modN = improve_EDDI_dictionary(ETo_7_day_modN, EDDI_dict_modN)

                
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_to_npy_file(dict EDDI_next_dict_modN, int model_NUM):
                
                    for idx_,i_val in enumerate(EDDI_next_dict_modN):
                      
                    #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict = EDDI_next_dict_modN[f'{i_val}'][0]
                        
                        for dic_init_and_eddi_val in EDDI_final_dict:
                            # print(dic_init_and_eddi_val)
                            #Open up the file and insert the value
                            init_day = list(dic_init_and_eddi_val.keys())[0]
                            
                            '''When appending using several scripts, sometimes the file will
                            load before another file has finished filling in the file and causing
                            a ValueError: cannot reshape array of size 15328 into shape (2,4,45,27,59)'''
                            
                            fileOut = f'EDDI_{init_day}.npy'
                            
                            try:
                                eddi_open = np.load(fileOut,allow_pickle=True)
                            except ValueError:
                                os.system('sleep 2')
                                eddi_open = np.load(fileOut,allow_pickle=True)
                            
                            lead_values = np.load(f'EDDI_{init_day}_julian_lead.npy',allow_pickle=True)
                            
                            #Sometimes we found some indexes that really don't matter 
                            try:
                                index_val=np.where(lead_values == int(i_val))[0][0]
                                eddi_open[0,index_val,i_Y,i_X] = list(dic_init_and_eddi_val.values())[0]
                                np.save(fileOut,eddi_open)
                            except IndexError:
                                np.save(fileOut,eddi_open)
                                pass
                    
                add_to_npy_file(EDDI_next_dict_modN,model_NUM)
                
                #release memory (didn't work)
                # del EDDI_next_dict_modN, EDDI_dict_modN, ETo_7_day_modN
                
    print(f'Completed date {_date} and saved into {home_dir}/EDDI_mod{model_NUM}.')
    #save the dates that were completed to not re-run
    if model_NUM == 3:
        os.system(f'echo Completed {_date} >> {script_dir}/EDDI_completed_npy.txt')
#%%
'''For some odd reason, this function will keep allocating new memory (even though
there is nothing visible in the function that I can remove). To get around this,
it appears that I can only do 9 total dates between all 3 of my processes before
memory gets too high (potential memory leak), so add a break'''

#_date=init_date_list[0]
'''Read EDDI_completed_npy.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/EDDI_completed_npy_{model_NAM1}.txt',dtype='str')
try:
    #first line contains a header, nothing with dates
    completed_dates = completed_dates[:,1]
except IndexError:
    completed_dates = ''
# completed_dates = pd.to_datetime(completed_dates[:],format='%Y-%m-%d')
#only work on dates that aren't completed

subset_completed_dates = [i[5:] for i in completed_dates]
new_directory = f'{home_dir}/EDDI_mod{model_NUM}/already_completed'
count=0
for _date in init_date_list[start_:end_]:    
    if _date[5:] not in subset_completed_dates:
        EDDI_function(start_, end_, model_NUM, init_date_list, _date,var)
        count+=1
        #save
        os.system('sleep 15') #sometimes saving seems to not fully capture 1 file
        os.system(f'cp {home_dir}/EDDI_mod{model_NUM}/*.npy {new_directory}/')
        if count == 2:
            print('Completed 2 dates.')
            break
        
    

# if __name__ == '__main__':
#     call_function(start_, end_, model_NUM, init_date_list, _date)
