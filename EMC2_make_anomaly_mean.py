#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly for RZSM.
    
Anomaly mean is calculated as the 7-day average of RZSM by 42-day window for each calendar day/lead date.

@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn

# dir1 = 'main_dir'
# start_ = int('start_init')
# end_ = start_ + int('init_step')
# model_NAM1 = 'model_name'
# var = 'variable_'  #for anomaly function

# Test for 1 step size and model
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'EMC'
var='RZSM'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
os.chdir(home_dir)

anom_dir = f'{home_dir}/anomaly'
new_dir_mean = f'{anom_dir}/mean_for_ACC' #for saving the mean value output
os.system(f'mkdir -p {new_dir_mean}')

save_MME_anomaly_mean = f'{home_dir}/anomaly/MME/mean_MME'
save_MME_anomaly = f'{home_dir}/anomaly/MME'

script_dir = f'{dir1}/Scripts'
save_complete_dir = f'{script_dir}/anomaly_mean_completed'
os.system(f'mkdir {save_complete_dir}')

#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

gridMET_dir = f'{dir1}/Data/gridMET'

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files

###Files
file_list = os.listdir()

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{home_dir}/{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var)    

'''Steps for anomaly calculation. 
1.) For each date:
    2.) Go thorugh all files and find the 1 week average for the same dates:
        2a.) Find the mean
        2b.) Save into its own file to later create the anomaly. This will likely
    be much quicker than initially
        
'''

            

#%%
'''Make empty text files to not have to re-run code'''
new_eto_anom = f'ETo_completed_mean_for_anomaly_{model_NAM1}.txt'

new_rzsm_anom = f'RZSM_completed_mean_for_anomaly_{model_NAM1}.txt'

#since this is a pretty slow process
for name in [new_eto_anom,new_rzsm_anom]:
    try:
        completed_dates = np.loadtxt(f'{script_dir}/{name}',dtype='str')         
    except OSError:
        os.system(f'touch {script_dir}/{name}')
        name_index = f'{name}'.split('_')[0]
        os.system(f'echo "Completed {name_index}" > {script_dir}/{name}')


#%%    
'''process SubX files and create EDDI values from Subx, soil moisture anomalies, and 
reference evapotranspiration anomalies'''
# def multiProcess_EDDI_SubX_TEST(_date):
def make_subX_anomaly(_date, var,HP_conus_mask,anomaly_spread, save_MME_anomaly_mean,save_MME_anomaly,init_date_list):
    # _date=init_date_list[0]
    # i_X,i_Y=10,10
    # anomaly_spread=42

        
    print(f'Calculating {var} mean on SubX for {_date} and saving as .nc4 in {new_dir_mean}.') 
    #Keep only files within a certain date (for faster processing of files instead of loading all of them into memory)
    def limit_file_list(_date):

        month_start = pd.to_datetime(_date).month
        
        dates_to_keep = []
        '''if within 3 months, keep files'''
        for file in sorted(glob(f'{home_dir}/{var}*.nc4')):
            test_1 = month_start - pd.to_datetime(file[-14:-4]).month
            test_2 = pd.to_datetime(file[-14:-4]).month - month_start
            
            if test_1 <=-9 or test_1>= 9:
                dates_to_keep.append(file)
            elif abs(test_1) <=3 or abs(test_2) <=3:
                dates_to_keep.append(file)
                
        return (dates_to_keep)
    
    dates_to_keep = limit_file_list(_date)
    
    #Open all files for faster processing
    if var == 'RZSM':
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested')
        #Use subx_all to create the full distribution of all data
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc4', combine='by_coords').persist()
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc4', concat_dim=['S'], combine='nested').persist() #Load all into memory (really big)
        # subx_all = xr.open_mfdataset(dates_to_keep, concat_dim=['S'], combine='nested', chunks={"S": 1},parallel='True').persist()

        subx_all = xr.open_mfdataset(dates_to_keep, concat_dim=['S'], combine='nested',parallel='True').persist() #load first
        
        #Get a file with the actual dates
        new_list = [i[-5:] for i in init_date_list]
        index_of_actual_date = new_list.index(_date[-5:])
        use_date = init_date_list[index_of_actual_date]
        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_{model_NAM1}_{use_date}.nc4', engine='netcdf4')
        
        
        #Test why there is an issue 
        os.getcwd()
        file_list = sorted(glob('soilw*'))
        
        for file in file_list:
            open_f = xr.open_dataset(file)
            open_f.shape()
        
    elif var == 'ETo':
        #Open all files for faster processing (for EDDI and ETo anomaly)
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested')
        subx_all = xr.open_mfdataset(dates_to_keep, concat_dim=['S'], combine='nested').persist()

        #Get a file with the actual dates
        new_list = [i[-5:] for i in init_date_list]
        index_of_actual_date = new_list.index(_date[-5:])
        use_date = init_date_list[index_of_actual_date]
        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_{use_date}.nc4')
        
    #Save outputs into this file
    out_template = xr.zeros_like(subx_out)

    #convert to a numpy array for testing
    test_arr = subx_all[list(subx_all.keys())[0]].to_numpy()
    
    test_arr.shape
    '''We can now select the values based on only the lead julian date values. Because
    files with nothing in the julian date for that day will have an np.nan
    
    Example: _date="2000-01-05"
    test_arr[:,0,1:42,10,10]    
    
    --first 4 julian days, [1,2,3,4] have no values. 5th julian day has a value.
    EMC has a 35 day forecast. Looks good
    
    array([       nan,        nan,        nan,        nan, 0.6368441 ,
           0.64035332, 0.64381248, 0.64739686, 0.65017569, 0.65021122,
           0.64838684, 0.64401394, 0.63913   , 0.6335783 , 0.65773058,
           0.62275946, 0.62674689, 0.63152933, 0.63385528, 0.63239986,
           0.63352191, 0.6384446 , 0.64490312, 0.65216124, 0.65724766,
           0.66248369, 0.66521734, 0.66108894, 0.64609063, 0.61645734,
           0.60894138, 0.64115506, 0.6683653 , 0.6811946 , 0.69890547,
           0.71548241, 0.72964418, 0.73898357, 0.74735188,        nan,
                  nan,        nan])
    
    '''
    
    
    
    for i_Y in range(subx_out[f'{var}'].shape[3]):
        for i_X in range(subx_out[f'{var}'].shape[4]):
            if _date == '2000-01-01':
                print(f'Working on lat {i_Y} and lon {i_X}')

            #This is for the mask of CONUS, don't do extra unneeded calculations
            if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                #%%
                def weekly_mean_of_file():
                    #append 7-day average from each ensemble member of each lead week from single file
                    all_mean_var_mod = {}
                    
                    for i in range(subx_out.model.shape[0]):
                        all_mean_var_mod[f'{i}'] = []

                    for mod in all_mean_var_mod.keys():
                        # print(mod)
                        for idx,julian_d in enumerate(subx_out.lead.values):
                            # print(idx)
                            try:
                                #Idx=0 is an instantaneous forecast. Not very skillful for ETo, but just as skill as week 1 with RZSM.
                                if idx<7:
                                    if idx==0 or idx==1:
                                        all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[f'{var}'].isel(lead=idx).isel(S=0, model=0, X=i_X, Y=i_Y).values)})
                                    else:
                                        #Take the averages that we can based on idx
                                        all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[f'{var}'].isel(lead=slice(idx-idx+1,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})
                                elif idx >= 7:
                                    #Subtract minus 6 because of indexing. First is leads 1-7, which equals 7 days
                                    all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[f'{var}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})
                            except IndexError:
                                pass                    
                                #Index error exists because there is no data at the end of the file
                                #That meets the current idx restrictions
                        
                    return(all_mean_var_mod)
                
                all_mean_var_mod=weekly_mean_of_file()
                 
                #%%
             
                def find_anomaly_with_weekly_mean(all_mean_var_mod,test_arr,anomaly_spread=42):
                    out_dict = {}
                    out_mean = {}
                    '''This will return all the data for each file for averages
                    because of day of year, we have to subset the +/- 42 days in 
                    a different way due to how the files day of year is setup'''
                   
                    for mod in all_mean_var_mod.keys():
                        # print(mod)
                        
                        out_mean[mod] = []
                        out_dict[mod] = []
                        for k, julian_d in enumerate((all_mean_var_mod[mod])):
                            julian_d = int(list(julian_d)[0])
                            '''Grab all days with 42 days of day of year, need
                            this weird approach because of slicing issues'''
                            if int(julian_d) <= anomaly_spread:
                                subtract_ = int(julian_d)-anomaly_spread
                                front_half = test_arr[:,int(mod),0:julian_d+anomaly_spread,i_Y,i_X]
                                front_half=np.concatenate((front_half),axis=0)
                                
                                back_half = test_arr[:,int(mod),subtract_:,i_Y,i_X] #take the last n values from the end of the list which is the latest julian dates
                                back_half=np.concatenate((back_half),axis=0)
                                
                                mean_value=bn.nanmean(np.concatenate((front_half,back_half)))

                            elif julian_d >= (366-anomaly_spread):
                                add_ = (366-julian_d)*-1
                                diff_ = anomaly_spread - add_
                                
                                #Example juilan day = 359
                                #sub_small1 = all days after julian day (only a few for this example (7 days))
                                sub_small1 = test_arr[:,int(mod),add_:,i_Y,i_X] 
                                sub_small1=np.concatenate((sub_small1))
                                #Now get the difference dates from the beginning of the time series
                                sub_small2 = test_arr[:,int(mod),0:diff_,i_Y,i_X] 
                                sub_small2=np.concatenate((sub_small2))

                                #Now get all the julian days before the current julian_day
                                second_lead = list(subx_all.lead.values).index(julian_d)
                                first_lead = list(subx_all.lead.values).index(julian_d-anomaly_spread)
                               
                                sub_small3 = test_arr[:,int(mod),first_lead:second_lead,i_Y,i_X] 
                                sub_small3=np.concatenate((sub_small3))
                                
                                mean_value=bn.nanmean(np.concatenate((sub_small1,sub_small2,sub_small3)))

                            else:
                                #Need to find a way to get the index value from subx_all lead dimension because numpy has no way to tell what the actual value is

                                first_lead = list(subx_all.lead.values).index(julian_d-anomaly_spread)
                                second_lead = list(subx_all.lead.values).index(julian_d+anomaly_spread)
                                #Seems easier just to get keep this loaded instead of changing to an array
                                sub_out = test_arr[:,int(mod),0:first_lead+second_lead,i_Y,i_X]
                                mean_value=bn.nanmean(sub_out)

                                '''for some reason, model 3 with ETo has some inf values
                                that are messing up anomaly calculation'''

                            out_mean[mod].append({str(julian_d):mean_value})
                                
                        
                            # julian_d_str = str(julian_d)
                            # out_dict[mod].append({list(julian_d.keys())[0]:list(julian_d.values())[0] - list(out_mean[mod][0].values())[0]})
                                
                                # 
                                
                                # try:
                                #     value_ = list(mean_var_modN[f'{julian_d}'][i].values())[0]
                                #     anomaly_list.append(value_ - mean_value)
                                # except AttributeError:
                                #     value_=mean_var_modN[f'{julian_d}'][i]
                                #     anomaly_list.append(value_ - mean_value)
                    return(out_mean)
                
                model_mean_vals =find_anomaly_with_weekly_mean(all_mean_var_mod,test_arr,anomaly_spread=42)
                    
                def add_mean_to_nc_file(model_mean_vals):
                    for mod in model_mean_vals.keys():
                        # model_mean_vals[mod]
                        
                        out_list = []
                        for k in model_mean_vals[mod]:
                            out_list.append((list(k.values())[0]))
                            
                        out_array = np.array(out_list)
    
                        #Add data to file
                        out_template[list(out_template.keys())[0]][:,int(mod),:,i_Y,i_X] = out_array
                    return(out_template)
                
                out_template = add_mean_to_nc_file(model_mean_vals)
                
                    
    #After all coords have been done, save the file
    
    fileOut = "{}/{}_mean_{}.nc4".format(new_dir_mean,var,_date)
    
    if var == 'RZSM':
        var_OUT = xr.Dataset(
            data_vars = dict(
                RZSM_anom = (['S', 'model','L','Y','X'], out_template[list(out_template.keys())[0]].values),
            ),
            coords = dict(
                X = out_template.X.values,
                Y = out_template.Y.values,
                L = np.arange(len(out_template.lead.values)),
                M = out_template.model.values,
                S = out_template.S.values,
            ),
            attrs = dict(
                Description = 'RZSM mean from 42 day window',)
        )  
    elif var == 'ETo':
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_anom = (['S', 'model','L','Y','X'], out_template[list(out_template.keys())[0]].values),
            ),
            coords = dict(
                X = out_template.X.values,
                Y = out_template.Y.values,
                L = np.arange(len(out_template.lead.values)),
                M = out_template.model.values,
                S = out_template.S.values,
            ),
            attrs = dict(
                Description = 'ETo mean from 42 day window',)
        )  
    
    var_OUT.to_netcdf(path = fileOut, mode ='w', engine='scipy')
    
    print(f'Completed date {_date} and saved into {new_dir_mean}.')
    
    #save the dates that were completed to not re-run
    os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_mean_for_anomaly_{model_NAM1}.txt')
    return(0)
#%%
#Run function

vars_to_process = ['RZSM','ETo']
for var in vars_to_process:
    # _date=init_date_list[0]
    '''Read {var}_completed_anomaly_nc_.txt file to not have to re-run extra code'''
    completed_dates = np.loadtxt(f'{script_dir}/{var}_completed_mean_for_anomaly_{model_NAM1}.txt',dtype='str')
    try:
        #first line contains a header, nothing with dates
        completed_dates = completed_dates[:,1]
    except IndexError:
        completed_dates = ''
    
    
    #only work on dates that aren't completed
    subset_completed_dates = [i for i in completed_dates]
    
    #TODO: Need to create a new object that stores all of the yearly dates within a single year.
    
    #Step 1, use init date list and get only the months and the dates
    init_date_list[0][5:]
    month_day=[i[5:] for i in init_date_list]
    #Next, remove any duplicates
    month_day = sorted([*set(month_day)])
    
    #Only process files that aren't in list
    for idx, md in enumerate(month_day):
        new_date = f'2000-{md}'
        # print(idx,md)
        if new_date not in subset_completed_dates:
            make_subX_anomaly(_date=new_date,var=var, HP_conus_mask = HP_conus_mask, anomaly_spread=42,\
                              save_MME_anomaly_mean=save_MME_anomaly_mean,save_MME_anomaly=save_MME_anomaly,init_date_list=init_date_list)
    

