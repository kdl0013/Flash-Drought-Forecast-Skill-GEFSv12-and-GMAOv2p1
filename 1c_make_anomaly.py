#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    EMC calculate anomaly for RZSM.
    
Anomaly is calculated as the 7-day average of RZSM by 7-day window for each lead time.
Grab all sample values from within a 42 day window to construct distribution.
Then subtract the mean from the value to find the anomaly.

EMC is already in volumetric soil moisture content m3/m3

This script only has to process 1 years worth of files because you need all files 
to create the distribution. While processing, it adds to all other files.

Each member realization was kept seperately from the others to look at how
each model is compared to each other. We also construct the multi-member emsemble
and save that into another file


Save lead files to specific format to only have weekly leads (otherwise the entire file
                                                              is mainly empty)
out_lead = np.array([0,1,2,3,4,5,6,3.4,3.5,3.6,4.6])


@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn

dir1 = 'main_dir'
start_ = int('start_init')
end_ = start_ + int('init_step')
model_NAM1 = 'model_name'
var = 'RZSM'  #for anomaly function

if var == "RZSM":
    varname='soilw_bgrnd'


# # Test for 1 step size and model
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
start_ = int('0')
end_ = start_ + int('40')
model_NAM1 = 'EMC'
var = 'soilw_bgrnd'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
os.chdir(home_dir)

rzsm_dir = f'{home_dir}' #don't have to change directory, no conversion done

anom_dir = f'{home_dir}/anomaly'
new_dir_mean = f'{anom_dir}/mean_for_ACC' #for saving the mean value output
os.system(f'mkdir -p {new_dir_mean}')

save_MME_anomaly_mean = f'{home_dir}/anomaly/MME/mean_MME'
save_MME_anomaly = f'{home_dir}/anomaly/MME'

script_dir = f'{dir1}/Scripts'

#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

merra_dir = f'{dir1}/Data/MERRA2_NASA_POWER'

#Mask for CONUS
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask files

###Files
file_list = os.listdir()

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    

def return_date_list(varname):
    date_list = []
    for file in sorted(glob(f'{rzsm_dir}/{varname}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(varname)    #get all dates of files

#Make a full time series (some models have missing dates)

'''Steps for anomaly calculation. 
1.) For each date:
    2.) Go thorugh all files and find the 1 week average for the same dates:
        2a.) Find the mean
        3.) Calculate anomaly
        4.) Append all subsequent year's files since we have the first day 
    (only need to process the first year)
        
'''

#%%    
'''process SubX files and create EDDI values from Subx, soil moisture anomalies, and 
reference evapotranspiration anomalies'''
# def multiProcess_EDDI_SubX_TEST(_date):
def make_subX_anomaly(_date, var,HP_conus_mask,anomaly_spread, save_MME_anomaly_mean,save_MME_anomaly):
    # _date=init_date_list[0]
    # i_X,i_Y=10,10
    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .nc4 in {anom_dir}.') 

    #Open all files for faster processing
    if var == 'RZSM':
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested')
        #Use subx_all to create the full distribution of all data
        subx_all = xr.open_mfdataset(f'{home_dir}/{varname}*.nc4', concat_dim=['S'], combine='nested').persist()
        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{varname}_{_date}.nc4')
        variable='SM_SubX_m3_m3'
        var_name = f'{variable}_value'

    elif var == 'ETo':
        #Open all files for faster processing (for EDDI and ETo anomaly)
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested')
        subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested').persist()

        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_{_date}.nc')
        var_name = var

    for i_Y in range(subx_out[f'{var_name}'].shape[3]):
        for i_X in range(subx_out[f'{var_name}'].shape[4]):
            if _date == '2000-01-05':
                print(f'Working on lat {i_Y} and lon {i_X}')

            '''Initially, I just wanted to use the SMERGE files as a mask, but it
            does unneccesary computations for a lot of grid cells
            
            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE. But just use the USDM Conus mask. Values 1-6 indicate a region with USDM mask.
            
            '''
            #This is for the mask of CONUS, don't do extra unneeded calculations
            if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                #%%
                def weekly_mean_by_lead_for_all_years():
                    #append 7-day average from each ensemble member of each lead week from single file
                    mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3 = {},{},{},{}

                    for idx,julian_d in enumerate(subx_out.lead.values):
                        try:
                            #Idx=0 is an instantaneous forecast. Not very skillful for ETo, but just as skill as week 1 with RZSM.
                            if idx == 0:
                                mean_var_mod0[f'{julian_d}']=[]
                                mean_var_mod0[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod1[f'{julian_d}']=[]
                                mean_var_mod1[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod2[f'{julian_d}']=[]
                                mean_var_mod2[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod3[f'{julian_d}']=[]
                                mean_var_mod3[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                            elif idx >= 7 :
                                #Subtract minus 6 because of indexing. First is leads 1-7, which equals 7 days
                                mean_var_mod0[f'{julian_d}']=[]
                                mean_var_mod0[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod1[f'{julian_d}']=[]
                                mean_var_mod1[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod2[f'{julian_d}']=[]
                                mean_var_mod2[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod3[f'{julian_d}']=[]
                                mean_var_mod3[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})
                            
                            ''''#if we want to look at weekly summation of ETo. But it's much harder to
                            #create an anomaly distribution if you apply a rolling mean to 
                            #a 45-day time series'''
                          
                        except IndexError:
                            pass                    
                            #Index error exists because there is no data at the end of the file
                            #That meets the current idx restrictions
                        
                    return(mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3)
                
                mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3=weekly_mean_by_lead_for_all_years()
                
                def collect_all_weekly_lead_averages_from_other_files(mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3):
                    '''7-day mean of RZSM or ETo by lead date for the first opened file:
                    Next we will append to each julian day value in the key in the dictionary with
                    the same julian day from all files'''
                    dates_to_keep = []
                    if var == 'RZSM':
                        '''grab yearly week lead files since we already have the distribution'''
                        for file in sorted(glob(f'{home_dir}/{var}_SubX*{_date[-5:]}.nc4')):
                                dates_to_keep.append(file)
                        day_s = -14
                        day_e = -4
                                
                    elif var == 'ETo':
                        for file in sorted(glob(f'{home_dir}/ETo*{_date[-5:]}.nc')):
                            #I have ETo_anomaly files also in the directory, don't include those
                            if 'anomaly' not in file:
                                dates_to_keep.append(file)
                            
                        day_s = -13
                        day_e = -3
                            
                    for file in dates_to_keep:
                        #Dont' re-open the same file
                        if file[day_s:day_e] != _date:
                            open_f = xr.open_dataset(file)
    
                            for idx,julian_d in enumerate(subx_out.lead.values):
                                #Only look at idx up to 39 because we need a full 7 days of data in order to calculate EDDI
                                #first day cannot have a weekly mean. Just take a single value (leave np.nanmean function so array comes out properly)
    
                                if idx == 0:
                                    mean_var_mod0[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod1[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod2[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod3[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
    
                            
                                elif idx >= 7 :
                                    try:
                                        mean_var_mod0[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod1[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod2[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod3[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-6,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
    
                                    #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass
                    return(mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3,dates_to_keep,day_s,day_e)
                
                mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3,dates_to_keep,day_s,day_e \
                    =collect_all_weekly_lead_averages_from_other_files(mean_var_mod0,mean_var_mod1,mean_var_mod2,mean_var_mod3)                    
                #%%
                def find_anomaly_with_weekly_mean(mean_var_modN,mod_number,anomaly_spread,day_s,day_e):
                    out_dict = {}
                    out_mean = {}
                    '''This will return all the data for each file for averages
                    because of day of year, we have to subset the +/- 42 days in 
                    a different way due to how the files day of year is setup'''
                   
                    for k, julian_d in enumerate(mean_var_modN):
                        # if out_lead[k] == 1 or out_lead[k] == 2 or out_lead[k] == 3 or out_lead[k] == 4 or k==5 or k==6 or k==7:
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

                        '''for some reason, model 3 with ETo has some inf values
                        that are messing up anomaly calculation'''
                        sub_out_vals = sub_out.values
                        sub_out_vals = sub_out_vals.flatten()
                        sub_out_vals[sub_out_vals > 1e06]  = np.nan

                        mean_value = np.nanmean(sub_out_vals)
                        
                        out_mean[f'{k}'] = mean_value
                        
                        #subtract mean
                        anomaly_list = []
                        julian_d = str(julian_d)
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                            try:
                                value_ = list(mean_var_modN[f'{julian_d}'][i].values())[0]
                                anomaly_list.append(value_ - mean_value)
                            except AttributeError:
                                value_=mean_var_modN[f'{julian_d}'][i]
                                anomaly_list.append(value_ - mean_value)

                        out_dict[f'{k}']=anomaly_list
                    #END of loop
                    out_mean_final = out_mean
                    
                    '''Find new means for weekly averages first
                    Averages are from all of the individual observations'''
                    out_mean_final['3.4'] = ((out_mean['3']+ out_mean['4'])/2).astype('float32')
                    out_mean_final['3.5'] = ((out_mean['3']+ out_mean['4'] + out_mean['5']) /3).astype('float32')
                    out_mean_final['3.6'] = ((out_mean['3']+ out_mean['4'] + out_mean['5']+ out_mean['6']) /4).astype('float32')
                    out_mean_final['4.6'] = ((out_mean['4']+ out_mean['5'] + out_mean['6']) /3).astype('float32')
                    
                    #get week lead first from Summation_ETo_modN, then subtract the new mean anomaly
                    
                    def create_multi_weekly_anom (out_dict,out_mean_final):
                        anom_list = []
                        # if week_num == '3&4':
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                            anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[3-1]][i].values())[0] + \
                                list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0])/2) - out_mean_final['3.4'])
                                
                            out_dict['3.4'] = anom_list 
                                

                        anom_list = []
                        # elif week_num == '3&4&5':
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                           anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[3-1]][i].values())[0] + \
                                list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0] + \
                                             list(mean_var_modN[list(mean_var_modN.keys())[5-1]][i].values())[0])/3) - out_mean_final['3.5'])
                           
                           out_dict['3.5'] = anom_list  
                        
                        anom_list = []
                        # elif week_num == '3&4&56':
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                           anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[3-1]][i].values())[0] + \
                                list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0] + \
                                             list(mean_var_modN[list(mean_var_modN.keys())[5-1]][i].values())[0] + \
                                                 list(mean_var_modN[list(mean_var_modN.keys())[6-1]][i].values())[0])/ 4) - out_mean_final['3.6'])
                           
                           out_dict['3.6'] = anom_list  
                        
                        anom_list = []
                        # elif week_num == '4&5&6':
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                           anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0] + \
                                list(mean_var_modN[list(mean_var_modN.keys())[5-1]][i].values())[0] + \
                                             list(mean_var_modN[list(mean_var_modN.keys())[6-1]][i].values())[0])/3) - out_mean_final['4.6'])
                           
                           out_dict['4.6'] = anom_list 
                           
                           
                           return(out_dict)
                           
                   #Call function
                    out_dict = create_multi_weekly_anom (out_dict=out_dict, out_mean_final=out_mean_final)

                    '''For some reason, the last weekly average lead doesn't work, just type manually'''
                    anom_list = []
                    # elif week_num == '4&5&6':
                    for i in range(len(mean_var_modN[f'{julian_d}'])):
                       anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0] + \
                            list(mean_var_modN[list(mean_var_modN.keys())[5-1]][i].values())[0] + \
                                         list(mean_var_modN[list(mean_var_modN.keys())[6-1]][i].values())[0])/3) - out_mean_final['4.6'])
                       
                       out_dict['4.6'] = anom_list 
                       
                    #End of find_anomaly_with_weekly_mean function
                    return (out_dict,out_mean_final)
                
                
                anom_mod0,mean_mod0 = find_anomaly_with_weekly_mean(mean_var_mod0,mod_number=0,anomaly_spread=anomaly_spread,day_s=day_s,day_e=day_e)
                anom_mod1,mean_mod1 = find_anomaly_with_weekly_mean(mean_var_mod1,mod_number=1,anomaly_spread=anomaly_spread,day_s=day_s,day_e=day_e)
                anom_mod2,mean_mod2 = find_anomaly_with_weekly_mean(mean_var_mod2,mod_number=2,anomaly_spread=anomaly_spread,day_s=day_s,day_e=day_e)
                anom_mod3,mean_mod3 = find_anomaly_with_weekly_mean(mean_var_mod3,mod_number=3,anomaly_spread=anomaly_spread,day_s=day_s,day_e=day_e)
#%%
                add_dates = [i[day_s:] for i in dates_to_keep]
                anom_mod0['init_dates'] = add_dates
                anom_mod1['init_dates'] = add_dates
                anom_mod2['init_dates'] = add_dates
                anom_mod3['init_dates'] = add_dates
                
                mean_avg_ensemble = {}
                for k,v in mean_mod0.items():
                    mean_avg_ensemble[f'{k}'] = (mean_mod0[f'{k}'] + mean_mod1[f'{k}'] + mean_mod2[f'{k}'] + mean_mod3[f'{k}'])/4
                
                mean_avg_ensemble['init_dates'] = add_dates #Save these values to the file for multi-model ensemble
                #%%
                #Create weekly leads all in one file
                '''We have created weekly leads for several weeks (see below)
                out_lead = np.array([0,1,2,3,4,5,6,3.4,3.5,3.6,4.6])
                
                Weeks: 1,2,3,4. Then we take averages of multiple weeks including:
                    weeks 3&4, 3&4&5, 3&4&5&6, 4&5&6.

                '''
                
                #Add all date files to it's own dictionary to quickly add to new files
                #This code will add each lead time to it's correct file (verified by my own eyes)
                def improve_anomaly_dictionary(anom_modN, dates_to_keep):
                    #first add a dictinary with the start date as the key
                    mod_out = {}
                    for idx in range(len(dates_to_keep)):
                        mod_out[anom_modN['init_dates'][idx]] = []
                        
                        all_values = []
                        count=0
                        for init_day,values in anom_modN.items():

                            if init_day=='init_dates':
                                break
                            else:
                                all_values.append(anom_modN[init_day][idx])
                                count+=1
                        mod_out[anom_modN['init_dates'][idx]] = all_values    

                    return(mod_out)
                
                out_mod0 = improve_anomaly_dictionary(anom_mod0, dates_to_keep)
                out_mod1 = improve_anomaly_dictionary(anom_mod1, dates_to_keep)
                out_mod2 = improve_anomaly_dictionary(anom_mod2, dates_to_keep)
                out_mod3 = improve_anomaly_dictionary(anom_mod3, dates_to_keep)

                
                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to all files.
                
                Since we are doing a multimodel ensemble, we will now take the average of the
                anomalies. (We already have the mean average of all models saved'''
                
                all_model_anomaly_mean = {}
                for idx,date in enumerate(out_mod0):
                    all_model_anomaly_mean[date]=list((np.array(out_mod0[date]) + \
                                                      np.array(out_mod1[date]) + \
                                                    np.array(out_mod2[date]) + \
                                                        np.array(out_mod3[date])) /4)

                def improve_MME_mean_dict(mean_avg_ensemble,dates_to_keep):
                    mod_out = {}
                    for idx in range(len(dates_to_keep)):
                        mod_out[mean_avg_ensemble['init_dates'][idx]] = []
                        
                        all_values = []
                        count=0
                        for init_day,values in mean_avg_ensemble.items():

                            if init_day=='init_dates':
                                break
                            else:
                                all_values.append(mean_avg_ensemble[init_day])
                                count+=1
                        mod_out[mean_avg_ensemble['init_dates'][idx]] = all_values   
                        
                    return(mod_out)
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
                
                add_MME_to_file = improve_MME_mean_dict(mean_avg_ensemble,dates_to_keep)
                '''Now that we have created new files, we can append each file with the data that was found'''
                
                def add_anomaly_to_nc_file(out_mod0,out_mod1,out_mod2,out_mod3):

                    for idx_,init_file_day in enumerate(out_mod0):
                        
                        #We have the file, now open it
                        #add values by julian day
                        if var == 'ETo':
                            fileOut = ("{}/{}_anomaly_{}4".format(anom_dir,var,init_file_day))
                        elif var == 'RZSM':
                            fileOut = ("{}/{}_anomaly_{}".format(anom_dir,var,init_file_day))
                            
                        var_n = f'{var}_anom'
                        file_open = xr.open_dataset(fileOut)
                        file_open.close()
                        
                        file_open[var_n][0,0,:,i_Y,i_X] = np.array(out_mod0[init_file_day])
                        file_open[var_n][0,1,:,i_Y,i_X] = np.array(out_mod1[init_file_day])
                        file_open[var_n][0,2,:,i_Y,i_X] = np.array(out_mod2[init_file_day])
                        file_open[var_n][0,3,:,i_Y,i_X] = np.array(out_mod3[init_file_day])

                        
                        file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                        file_open.close()

                add_anomaly_to_nc_file(out_mod0, out_mod1,out_mod2,out_mod3)
                
                def add_mean_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,dates_to_keep):

                    for idx_,lead in enumerate(mean_mod0):
                        #We have the mean value for each year/lead time/model
                        #Add to only 1 file because its the mean of all years
                                                
                        if var == 'RZSM':
                            day_s = -14
                            fileOut = "{}/{}_mean_{}".format(new_dir_mean,var,dates_to_keep[0][day_s:])

                        elif var == 'ETo':
                            day_s = -13
                            fileOut = "{}/{}_mean_{}4".format(new_dir_mean,var,dates_to_keep[0][day_s:])

                        file_open = xr.open_dataset(fileOut)
                        file_open.close()
                        
                        var_n = f'{var}_mean'
                        
                        #Add data to netcdf file
                        file_open[var_n][0,0,idx_,i_Y,i_X] = mean_mod0[f'{lead}']
                        file_open[var_n][0,1,idx_,i_Y,i_X] = mean_mod1[f'{lead}']
                        file_open[var_n][0,2,idx_,i_Y,i_X] = mean_mod2[f'{lead}']
                        file_open[var_n][0,3,idx_,i_Y,i_X] = mean_mod3[f'{lead}']
                        
                        file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                        file_open.close()

                        
                add_mean_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,dates_to_keep)
                
                #ADD MME to files
                def add_anomaly_to_nc_file_MME(all_model_anomaly_mean):

                    for idx_,init_file_day in enumerate(all_model_anomaly_mean):
                        
                        #We have the file, now open it
                        #add values by julian day
                        if var == 'ETo':
                            fileOut = ("{}/{}_MME_anomaly_{}4".format(save_MME_anomaly,var,init_file_day))
                        elif var == 'RZSM':
                            fileOut = ("{}/{}_MME_anomaly_{}".format(save_MME_anomaly,var,init_file_day))
                            
                        var_n = f'{var}_MME_anom'
                        file_open = xr.open_dataset(fileOut)
                        file_open.close()
                        #Now add to file by lead week (or average)
                        file_open[var_n][0,:,i_Y,i_X] = np.array(all_model_anomaly_mean[init_file_day])
                        file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                        file_open.close()
                        
                     
                add_anomaly_to_nc_file_MME(all_model_anomaly_mean)
                
                
                
                def add_mean_to_nc_file_MME(add_MME_to_file):
                    
                    #Only need to add the first date to the file (since we only need 1 year's worth of data as the mean)
                    for idx_,init_file_day in enumerate(add_MME_to_file):
                        #We have the mean value for each year/lead time/model
                        #Add to only 1 file because its the mean of all years
                                                
                        if var == 'RZSM':
                            fileOut = "{}/{}_MME_mean_{}.nc4".format(save_MME_anomaly_mean,var,_date)

                        elif var == 'ETo':
                            fileOut = "{}/{}_MME_mean_{}.nc4".format(save_MME_anomaly_mean,var,_date)

                        file_open = xr.open_dataset(fileOut)
                        file_open.close()
                        
                        var_n = f'{var}_MME_mean'
                        
                        lead_values = add_MME_to_file[f'{_date}.nc4']       
                        file_open[var_n][0,:,i_Y,i_X] = np.array(lead_values)
                        
                        file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                        file_open.close()
                        
                add_mean_to_nc_file_MME(add_MME_to_file)
#%%
                
    print(f'Completed date {_date} and saved into {home_dir}.')
    
    #save the dates that were completed to not re-run
    os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_anomaly_nc_{model_NAM1}.txt')

    return()
#%%



# _date=init_date_list[0]
'''Read {var}_completed_anomaly_nc_.txt file to not have to re-run extra code'''
completed_dates = np.loadtxt(f'{script_dir}/{var}_completed_anomaly_nc_{model_NAM1}.txt',dtype='str')
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
        make_subX_anomaly(_date=_date,var=var, HP_conus_mask = HP_conus_mask, anomaly_spread=42,\
                          save_MME_anomaly_mean=save_MME_anomaly_mean,save_MME_anomaly=save_MME_anomaly)



# count=0
# for _date in init_date_list[start_:end_]:    
#     if start_ == 50 and _date == '2000-01-10':
#         break
#     else:
#         RZSM_anomaly(start_,end_,init_date_list, _date,var)


