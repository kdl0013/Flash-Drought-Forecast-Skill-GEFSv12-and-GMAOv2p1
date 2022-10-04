#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate anomaly for ETo (no RZSM).
    
Anomaly is calculated as the 7-day average of ETo by 42-day window for each lead time.



@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn
from multiprocessing import Pool

# dir1 = 'main_dir'
# start_ = int('start_init')
# end_ = start_ + int('init_step')
# model_NAM1 = 'model_name'
# var = 'variable_'  #for anomaly function

# Test for 1 step size and model
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'model_name'
var = 'ETo'

if model_NAM1 == "GMAO":
    n_processes = 3
else:
    n_processes = 2

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'

anom_dir = f'{home_dir}/anomaly'
new_dir_mean = f'{anom_dir}/mean_for_ACC' #for saving the mean value output
os.system(f'mkdir -p {new_dir_mean}')

save_MME_anomaly_mean = f'{home_dir}/anomaly/MME/mean_MME'
save_MME_anomaly = f'{home_dir}/anomaly/MME'

script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)
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
'''process SubX files and create EDDI values from Subx, soil moisture anomalies, and 
reference evapotranspiration anomalies'''
# def multiProcess_EDDI_SubX_TEST(_date):
def make_subX_anomaly(_date):

    # _date=init_date_list[0]
    # i_X,i_Y=10,10
    # anomaly_spread=42
    #TODO: Change to penman or priestley for name
    name_='Penman'  #or Penman
    fileOut = "{}/{}_mean_{}_{}.nc4".format(new_dir_mean,var,name_,_date)

    try:
        xr.open_dataset(f'{fileOut}',engine='netcdf4')
    except FileNotFoundError:
            
    
        print(f'Calculating {var} mean on SubX for {_date} and saving as .nc4 in {new_dir_mean}.') 
        #Keep only files within a certain date (for faster processing of files instead of loading all of them into memory)
        def limit_file_list(_date):
    
            month_start = pd.to_datetime(_date).month
            
            dates_to_keep = []
            
            '''if within 3 months, keep files'''
            for file in sorted(glob(f'{home_dir}/{var}*{name_}*.nc4')):
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
    
            subx_all = xr.open_mfdataset(dates_to_keep, concat_dim=['S'], combine='nested',parallel='True') #load first
            
            #Get a file with the actual dates
            new_list = [i[-5:] for i in init_date_list]
            index_of_actual_date = new_list.index(_date[-5:])
            use_date = init_date_list[index_of_actual_date]
            #Process indivdual files for output
            subx_out = xr.open_dataset(f'{home_dir}/{var}_{name_}_{model_NAM1}_{use_date}.nc4', engine='netcdf4')
                        
        elif var == 'ETo':
            #Open all files for faster processing (for EDDI and ETo anomaly)
            # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested')
            subx_all = xr.open_mfdataset(dates_to_keep, concat_dim=['S'], combine='nested',parallel='True')
            # tt=subx_all.ETo[:,0,:,10,10].values
            #Get a file with the actual dates
            new_list = [i[-5:] for i in init_date_list]
            index_of_actual_date = new_list.index(_date[-5:])
            use_date = init_date_list[index_of_actual_date]
            #Process indivdual files for output
            subx_out = xr.open_dataset(f'{home_dir}/{var}_{name_}_{model_NAM1}_{use_date}.nc4')
            
            '''Fixed code with if statement if there are different number of models'''
            if model_NAM1 == 'ESRL':
                if subx_all.M.shape[0] != subx_out.M.shape[0]:
                    subx_all = subx_all.sel(M=slice(1,5))
            
            
        '''Test why date with ESRL 2000-03-01 produces 5 models when there are only 4'''
        # os.chdir('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/ESRL')
        # out_ = []
        # for f in sorted(glob('ETo_Prie*')):
        #     open_f = xr.open_dataset(f)
        #     out_.append(open_f.M.shape[0])            
        # np.unique(out_)
        
        # out_2 = []
        # for f in dates_to_keep:
        #     open_f = xr.open_dataset(f)
        #     out_2.append(open_f.M.shape[0])            
        # np.unique(out_2)
        # '''All files have 4 models, not sure why a 5th model is introduced and messing up the code.
        # Test to see what models have values'''
        # #%%
        # for moo in range(subx_all.M.shape[0]):
        #     print(subx_all.ETo.isel(M=moo,Y=10,X=10,L=slice(4,10)).values)
        # '''After testing, the 0th model has no values. Fix above code with if statement'''
        
            
            
        #Save outputs into this file
        out_template = xr.zeros_like(subx_out)
    
        #convert to a numpy array for testing
        test_arr = subx_all[list(subx_all.keys())[0]].to_numpy()
        # del test_arr        
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
        
        
        
        for i_Y in range(subx_out.Y.shape[0]):
            for i_X in range(subx_out.X.shape[0]):
                if _date == '2000-01-01':
                    print(f'Working on lat {i_Y} and lon {i_X}')
    
                #This is for the mask of CONUS, don't do extra unneeded calculations
                if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                    #%%
                    def weekly_mean_of_file():
                        #append 7-day average from each ensemble member of each lead week from single file
                        all_mean_var_mod = {}
                        
                        for i in range(subx_out.M.shape[0]):
                            all_mean_var_mod[f'{i}'] = []
    
                        for mod in all_mean_var_mod.keys():
                            # print(mod)
                            for idx,julian_d in enumerate(subx_out.L.values):
                                # print(idx)
                                try:
                                    #Idx=0 is an instantaneous forecast. Not very skillful for ETo, but just as skill as week 1 with RZSM.
                                    if idx<7:
                                        if idx==0 or idx==1:
                                            all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[list(subx_out.keys())[0]].isel(L=idx).isel(S=0, M=int(mod), X=i_X, Y=i_Y).values)})
                                        else:
                                            #Take the averages that we can based on idx
                                            all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[list(subx_out.keys())[0]].isel(L=slice(idx-idx+1,idx)).isel(S=0, M=int(mod), X=i_X, Y=i_Y).values)})
                                    elif idx >= 7:
                                        #Subtract minus 6 because of indexing. First is leads 1-7, which equals 7 days
                                        all_mean_var_mod[mod].append({f'{julian_d}':bn.nanmean(subx_out[list(subx_out.keys())[0]].isel(L=slice(idx-6,idx)).isel(S=0, M=int(mod), X=i_X, Y=i_Y).values)})
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
                                # if julian_d == 366:
                                #     break
                                '''Grab all days with 42 days of day of year, need
                                this weird approach because of slicing issues'''
                                #If doy <= 42
                                if int(julian_d) < anomaly_spread:
                                    subtract_ = int(julian_d)-anomaly_spread #what to grabl from back of dataset
                                    front_half = test_arr[:,int(mod),0:julian_d+anomaly_spread,i_Y,i_X] # doy 1 through front of dataset
                                    front_half=np.concatenate((front_half),axis=0) #flatten
                                    
                                    back_half = test_arr[:,int(mod),subtract_:,i_Y,i_X] #take the last n values from the end of the list which is the latest julian dates
                                    back_half=np.concatenate((back_half),axis=0)
                                    
                                    mean_value=bn.nanmean(np.concatenate((front_half,back_half)))
                                #If doy >= 324
                                elif julian_d > (366-anomaly_spread) and julian_d != 366:
    
                                    # break
                                    add_ = (366-julian_d)*-1 #going to get from the front half of dataset
                                    diff_ = anomaly_spread - add_ #total other days missing
                                    
                                    #Example juilan day = 359
                                    #sub_small1 = all days after julian day (only a few for this example (7 days))
                                    sub_small1 = test_arr[:,int(mod),add_:,i_Y,i_X] 
                                    sub_small1=np.concatenate((sub_small1))
                                    #Now get the difference dates from the beginning of the time series
                                    sub_small2 = test_arr[:,int(mod),0:diff_,i_Y,i_X] 
                                    sub_small2=np.concatenate((sub_small2))
    
                                    #Now get all the julian days before the current julian_day
                                    second_lead = list(subx_all.L.values).index(julian_d)
                                    first_lead = list(subx_all.L.values).index(julian_d-anomaly_spread)
                                   
                                    sub_small3 = test_arr[:,int(mod),first_lead:second_lead,i_Y,i_X] 
                                    sub_small3=np.concatenate((sub_small3))
                                    
                                    mean_value=bn.nanmean(np.concatenate((sub_small1,sub_small2,sub_small3)))
    
                                else:
                                    if julian_d == 366:
                                        sub_small1 = test_arr[:,int(mod),0:anomaly_spread,i_Y,i_X]
                                        sub_small1=np.concatenate((sub_small1))
    
                                        sub_small2 = test_arr[:,int(mod),-anomaly_spread:,i_Y,i_X]
                                        sub_small2=np.concatenate((sub_small2))
                                        mean_value=bn.nanmean(np.concatenate((sub_small1,sub_small2)))
                                    else:
                                        #Need to find a way to get the index value from subx_all lead dimension because numpy has no way to tell what the actual value is
                                        if julian_d ==anomaly_spread:
                                            #fix indexing issues 
                                            first_lead = list(subx_all.L.values).index(julian_d+1-anomaly_spread) 
                                            second_lead = list(subx_all.L.values).index(julian_d+1+anomaly_spread)
                                        else:
                                            first_lead = list(subx_all.L.values).index(julian_d-anomaly_spread) 
                                            second_lead = list(subx_all.L.values).index(julian_d+anomaly_spread)
    
                                        #Seems easier just to get keep this loaded instead of changing to an array
                                        sub_out = test_arr[:,int(mod),first_lead:second_lead,i_Y,i_X]
                                        mean_value=bn.nanmean(sub_out)
    
                                out_mean[mod].append({str(julian_d):mean_value})
    
                        return(out_mean)
                    
                    model_mean_vals =find_anomaly_with_weekly_mean(all_mean_var_mod,test_arr,anomaly_spread=42)
                    model_mean_vals
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
    
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_mean = (['S', 'M','L','Y','X'], out_template[list(out_template.keys())[0]].values),
            ),
            coords = dict(
                X = out_template.X.values,
                Y = out_template.Y.values,
                L = np.arange(len(out_template.L.values)),
                M = out_template.M.values,
                S = out_template.S.values,
            ),
            attrs = dict(
                Description = f'ETo {name_} mean from 42 day window',)
        )  
        
        var_OUT.to_netcdf(path = fileOut, mode ='w', engine='scipy')
        
        print(f'Completed date {_date} and saved into {new_dir_mean}.')
        
    # #save the dates that were completed to not re-run
    # os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_mean_for_anomaly_{model_NAM1}.txt')
    return(0)
#%%
#Step 1, use init date list and get only the months and the dates for a single year
#as the mean file. Then we can subtract to make the anomaly later
init_date_list[0][5:]
month_day=[i[5:] for i in init_date_list]
#Next, remove any duplicates
month_day = sorted([*set(month_day)])
month_day=[f'2000-{i}' for i in month_day]

if __name__ == '__main__':
    p = Pool(n_processes)
    p.map(make_subX_anomaly,month_day)

#Run a loop through only 1 year of files
#Only process files that aren't in list
# for idx, md in enumerate(month_day):
#     new_date = f'2000-{md}'
#     # print(idx,md)
#     # if new_date not in subset_completed_dates:
#     make_subX_anomaly(_date=new_date,var=var, HP_conus_mask = HP_conus_mask, anomaly_spread=42,\
#                       save_MME_anomaly_mean=save_MME_anomaly_mean,save_MME_anomaly=save_MME_anomaly,init_date_list=init_date_list)


