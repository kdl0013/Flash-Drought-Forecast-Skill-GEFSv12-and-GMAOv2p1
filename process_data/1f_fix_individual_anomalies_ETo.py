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

import warnings
warnings.filterwarnings("ignore")

dir1 = 'main_dir'
model_NAM1 = 'model_name'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
script_dir = f'{dir1}/Scripts'

os.chdir(home_dir)

gridMET_dir = f'{dir1}/Data/gridMET'
smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask file
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


'''Now we want to open up each ETo- and RZSM anomaly file (also EDDI), to see
which ones have missing data'''
var = 'RZSM'
file_list = f'{var}_anomaly_*.nc4'
missing_anomaly_var = [] #subx missing anomaly data files
missing_original_var = []

# test_file = f'{home_dir}/test/ETo_anomaly_2000-09-17.nc'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    # open_f = xr.open_dataset(test_file)

    open_f.close()
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    
    #anomaly. I inspected a file that was filled correctly and it had 528 missing value
    #If file is completely empty, the file will have all zeros and a unique value of only 1
    if len(np.unique(open_f[f'{var}_anom'][:,:,7,:,:])) <= 1 or \
        len(np.unique(open_f[f'{var}_anom'][:,:,14,:,:])) <= 1 or \
            len(np.unique(open_f[f'{var}_anom'][:,:,21,:,:])) <= 1 or \
                len(np.unique(open_f[f'{var}_anom'][:,:,28,:,:])) <= 1 or \
                    len(np.unique(open_f[f'{var}_anom'][:,:,35,:,:])) <= 1 or \
                        len(np.unique(open_f[f'{var}_anom'][:,:,42,:,:])) <= 1:
        missing_anomaly_var.append(file)    

    '''I manually calculated ETo from SubX variables, look if I am missing any
    data from those files'''

    if var == 'ETo':
        open_ETo = xr.open_dataset(f'{var}_{file[-13:-3]}.nc4') 
        if np.count_nonzero(np.isnan(open_ETo.ETo[0,:,:,:,:].values)) == 286740:
            missing_original_var.append(file)
    elif var == 'RZSM':
        open_m3_m3 = xr.open_dataset(f'SM_converted_m3_m3/SM_SubX_m3_m3_{file[-14:-4]}.nc4') 
        if np.count_nonzero(np.isnan(open_m3_m3.SM_SubX_m3_m3_value[0,:,:,:,:].values)) == 286740:
            missing_original_var.append(file)


missing_anomaly_var[0][-14:-4]
missing_dates = [i[-14:-4] for i in missing_anomaly_var]

print(f'Missing data in {len(missing_dates)}. Issue is likely leap year index issues. Fixing now if needed.')
print(f'Missing {var} SubX data (no anomaly) in {len(missing_original_var)} files.')

print('Copying all files into new directory to process seperately for certain leap years because of julian day issues')
new_dir = 'ETo_anomaly_leap_year'

os.system(f'mkdir {new_dir}')
for file in sorted(glob(file_list)):
    os.system(f"cp {file.split('/')[-1]} {new_dir}/")

#Move into new directory to process files
os.chdir(f'{home_dir}/{new_dir}')
#%%    
'''process SubX files that were messed up. I think I have finally figured out all
of the issues that I was having with some of the files'''
# def multiProcess_EDDI_SubX_TEST(_date):
# _date=missing_dates[0]
def anomaly_fix(_date, var, HP_conus_mask):
    #_date=init_date_list[0]  
  
    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .nc in {home_dir}') 

    #get julian day, timestamps, and the datetime
    def date_file_info( _date, variable):
        
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
    
    subx,file_timestamp_list,file_julian_list,file_datetime,var_name = date_file_info(_date=_date,variable='ETo')
     
    #Now convert to julian date and append coordinates
    subx2 = subx.assign_coords(lead = file_julian_list)
    
    # def moving_average(a, n=7) :
    #     ret = np.cumsum(a, dtype=float)
    #     ret[n:] = ret[n:] - ret[:-n]
    #     return ret[n - 1:] / n
    
    #Test on small subset
    #for i_Y in range(1):
    for i_Y in range(subx2[f'{var_name}'].shape[3]):
        #for i_X in range(1):
         for i_X in range(subx2[f'{var_name}'].shape[4]):

            #only work on grid cells within CONUS
            if HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7):

                def dict1_subx2():
                    
                    #append 7-day average from all files to a new dictionary
                    summation_ETo_mod0 = {}
                    summation_ETo_mod1 = {}
                    summation_ETo_mod2 = {}
                    summation_ETo_mod3 = {}

                    
                    for idx,julian_d in enumerate(file_julian_list):

                        if idx % 7 == 0 and idx !=0:
                            try:
                                summation_ETo_mod0[f'{julian_d}']=[]
                                summation_ETo_mod0[f'{julian_d}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod1[f'{julian_d}']=[]
                                summation_ETo_mod1[f'{julian_d}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod2[f'{julian_d}']=[]
                                summation_ETo_mod2[f'{julian_d}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                summation_ETo_mod3[f'{julian_d}']=[]
                                summation_ETo_mod3[f'{julian_d}'].append({f'{_date}':np.nanmean(subx2[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})
                            except IndexError:
                                pass

                    '''7-day mean of EDDI by index:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''
                    
                    dates_to_keep = []
                    '''if within 3 months, keep files'''
                    len_text=23
                    #ETo_anomaly is also in this directory, we don't want those
                    for file in sorted(glob(f'{home_dir}/ETo*{_date[-5:]}.nc4')):
                        if len(file.split('/')[-1]) == len_text:
                            dates_to_keep.append(file)
                                        
                    # for idx,file in dates_to_keep:
                    #     if dates_to_keep[-1][-13:-3] == _date:
                    #         pass
                        #Dont' re-open the same file
                        if file[-14:-4] != _date:
                            #Open up ETo file
                            open_f = xr.open_dataset(file)
                            Et_ref_open_f = open_f.assign_coords(lead = file_julian_list)
                            
                            for idx,val in enumerate(file_julian_list):
                                
                                if idx % 7 == 0 and idx != 0:
                                    try:
                                        summation_ETo_mod0[f'{val}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod1[f'{val}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod2[f'{val}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        summation_ETo_mod3[f'{val}'].append({f'{file[-14:-4]}':np.nanmean(Et_ref_open_f.ETo.isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                                        
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
                def compute_anomaly(ETo_7_day_average_modN):
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
                    
                ETo_dict_mod0 = compute_anomaly(ETo_7_day_average_mod0)
                ETo_dict_mod1 = compute_anomaly(ETo_7_day_average_mod1)
                ETo_dict_mod2 = compute_anomaly(ETo_7_day_average_mod2)
                ETo_dict_mod3 = compute_anomaly(ETo_7_day_average_mod3)

                def improve_RZSM_dictionary(ETo_7_day_average_modN,  ETo_dict_modN):

                    
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
                
                ETo_next_dict_mod0 = improve_RZSM_dictionary(ETo_7_day_average_mod0, ETo_dict_mod0)
                ETo_next_dict_mod1 = improve_RZSM_dictionary(ETo_7_day_average_mod1, ETo_dict_mod1)
                ETo_next_dict_mod2 = improve_RZSM_dictionary(ETo_7_day_average_mod2, ETo_dict_mod2)
                ETo_next_dict_mod3 = improve_RZSM_dictionary(ETo_7_day_average_mod3, ETo_dict_mod3)
                
                def add_to_nc_file( ETo_next_dict_mod0, ETo_next_dict_mod1, ETo_next_dict_mod2, ETo_next_dict_mod3):

                    for idx_,i_val in enumerate(ETo_next_dict_mod0):
                      
                    #for some reason it's a list in a list, this fixes that, loop through julian day
                        EDDI_final_dict0 = ETo_next_dict_mod0[f'{i_val}'][0]
                        EDDI_final_dict1 = ETo_next_dict_mod1[f'{i_val}'][0]
                        EDDI_final_dict2 = ETo_next_dict_mod2[f'{i_val}'][0]
                        EDDI_final_dict3 = ETo_next_dict_mod3[f'{i_val}'][0]

                        for idx,dic_init_and_eddi_val in enumerate(EDDI_final_dict0):
                            #Only work on dates that are missing
                            if list(dic_init_and_eddi_val.items())[0][0] in missing_dates:
    
                                # print(dic_init_and_eddi_val)
                                #Open up the file and insert the value
                                init_day = list(dic_init_and_eddi_val.keys())[0]
                                
                                '''When appending using several scripts, sometimes the file will
                                load before another file has finished filling in the file and causing
                                a ValueError: cannot reshape array of size 15328 into shape (2,4,45,27,59)'''
                                var2 = 'ETo'
                                lead_values = np.load(f'{home_dir}/julian_lead_{init_day}.npy',allow_pickle=True)
                                
                                #add values by julian day
                                # fileOut = f'{var}_anomaly_{init_day}.npy'

                                fileOut = ("{}_anomaly_{}.nc".format(var2,init_day))
                                file_open = xr.open_dataset(fileOut)
                                file_open.close()
                                
                                index_val=np.where(lead_values == int(i_val))[0][0]
                                
                                
                                #Add data to netcdf file
                                file_open.ETo_anom[0,0,index_val,i_Y,i_X] = list(EDDI_final_dict0[idx].values())[0]
                                file_open.ETo_anom[0,1,index_val,i_Y,i_X] = list(EDDI_final_dict1[idx].values())[0]
                                file_open.ETo_anom[0,2,index_val,i_Y,i_X] = list(EDDI_final_dict2[idx].values())[0]
                                file_open.ETo_anom[0,3,index_val,i_Y,i_X] = list(EDDI_final_dict3[idx].values())[0]
       
                                file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                                file_open.close()
                              
                    
                add_to_nc_file(ETo_next_dict_mod0,ETo_next_dict_mod1,ETo_next_dict_mod2,ETo_next_dict_mod3)
             
    return()
#END FUNCTION
#%%
#Call function


#Run single date test
# _date = '1999-03-01'
# var
if len(missing_dates) == 0:
    print('There are no missing datafiles for {var}_anomaly from SubX')
else:
    '''Because we only need to process 1 year's worth of data, keep up with a list
        so that it doesn not re-run extra computation when it does not need to'''
    completed_dates = []  
    for _date in missing_dates:
        if _date[-5:] not in completed_dates:
            anomaly_fix(_date,var, HP_conus_mask)
            completed_dates.append(_date[-5:]) #don't re-iterate over completed files, it loops over all years
  
        
#Now move the missing date files back to home_dir
for _date in missing_dates:
    os.system(f"cp {new_dir}/ETo_anomaly_{_date}.nc {home_dir}/")
    
    
#%% Check again and see if now we have no files that are empty
os.chdir(home_dir)

var = 'ETo'
file_list = f'{var}_anomaly_*.nc'
missing_anomaly= [] #subx missing anomaly data files

# test_file = f'{home_dir}/test/ETo_anomaly_2000-09-17.nc'


for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    # open_f = xr.open_dataset(test_file)

    open_f.close()
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    #anomaly. I inspected a file that was filled correctly and it had 528 missing value
    #If file is completely empty, the file will have all zeros and a unique value of only 1
    if len(np.unique(open_f.ETo_anom[:,:,7,:,:])) <= 1 or \
        len(np.unique(open_f.ETo_anom[:,:,14,:,:])) <= 1 or \
            len(np.unique(open_f.ETo_anom[:,:,21,:,:])) <= 1 or \
                len(np.unique(open_f.ETo_anom[:,:,28,:,:])) <= 1 or \
                    len(np.unique(open_f.ETo_anom[:,:,35,:,:])) <= 1 or \
                        len(np.unique(open_f.ETo_anom[:,:,42,:,:])) <= 1:
        missing_anomaly.append(file)    


missing_anomaly[0][-13:-3]
missing_dates = [i[-13:-3] for i in missing_anomaly]

print(f'Missing data in {len(missing_anomaly)}. There should be no issues now.')
