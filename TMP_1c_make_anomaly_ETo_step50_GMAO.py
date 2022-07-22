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

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
start_ = int('50')
end_ = start_ + int('25')
model_NAM1 = 'GMAO'
var = 'ETo'  #for anomaly function

# # Test for 1 step size and model
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# start_ = int('0')
# end_ = start_ + int('40')
# model_NAM1 = 'GMAO'
# var = 'ETo'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/{model_NAM1}'
rzsm_dir = f'{home_dir}/SM_converted_m3_m3'

anom_dir = f'{home_dir}/anomaly'
new_dir_mean = f'{anom_dir}/mean_for_ACC' #for saving the mean value output
os.system(f'mkdir {home_dir}/{new_dir_mean}')


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

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
    
var_list='SM'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{rzsm_dir}/{var_list}*.nc4')):
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
        
'''

#%%    
'''process SubX files and create EDDI values from Subx, soil moisture anomalies, and 
reference evapotranspiration anomalies'''
# def multiProcess_EDDI_SubX_TEST(_date):
def make_subX_anomaly(_date, var,HP_conus_mask,anomaly_spread):
    # _date=init_date_list[0]
    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .nc4 in {anom_dir}.') 
    #Used for eliminating iterating over grid cells that don't matter
    smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
    smerge_file_julian = smerge_file.copy()
  
    #Open all files for faster processing
    if var == 'RZSM':
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested')
        subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested').persist()
        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_SubX_m3_m3_{_date}.nc4')
        variable='SM_SubX_m3_m3'
        var_name = f'{variable}_value'

    elif var == 'ETo':
        #Open all files for faster processing (for EDDI and ETo anomaly)
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested')
        subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested').persist()
        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_{_date}.nc')
        var_name = var

    #Convert to julian day for processing
    smerge_day = pd.to_datetime(smerge_file.CCI_ano.time.values)
    smerge_julian = [i.timetuple().tm_yday for i in smerge_day]
    
    smerge_file_julian= smerge_file_julian.assign_coords({'time': smerge_julian})
    
    for i_Y in range(subx_out[f'{var_name}'].shape[3]):
        for i_X in range(subx_out[f'{var_name}'].shape[4]):
            if _date == '1999-01-10':
                print(f'Working on lat {i_Y} and lon {i_X}')

            '''Initially, I just wanted to use the SMERGE files as a mask, but it
            does unneccesary computations for a lot of grid cells
            
            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE
            
            if var == 'RZSM':
                true_or_false = (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1)
            elif var == 'ETo':
                true_or_false = (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            '''
            #This is for the mask of CONUS, don't do extra unneeded calculations
            if (HP_conus_mask.High_Plains[0,i_Y,i_X].values in np.arange(1,7)):
                
                def dict1_subx2():
                    
                    #append 7-day average from all files to a new dictionary
                    mean_var_mod0 = {}
                    mean_var_mod1 = {}
                    mean_var_mod2 = {}
                    mean_var_mod3 = {}
                    
                    for idx,julian_d in enumerate(subx_out.lead.values):
                        #You must julian_d + week_lead because with RZSM you need 7-day looking backwards into the past. Must have 7 values.
                        #Choose just model = 0 because we just need to know if there are 7-days total in any model
                        # print(idx)

                        try:
                            
                            if idx == 0:
                                mean_var_mod0[f'{julian_d}']=[]
                                mean_var_mod0[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod1[f'{julian_d}']=[]
                                mean_var_mod1[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod2[f'{julian_d}']=[]
                                mean_var_mod2[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod3[f'{julian_d}']=[]
                                mean_var_mod3[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=idx).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   
                            
                            
                            if idx % 7 == 0 and idx !=0:
                                mean_var_mod0[f'{julian_d}']=[]
                                mean_var_mod0[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod1[f'{julian_d}']=[]
                                mean_var_mod1[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod2[f'{julian_d}']=[]
                                mean_var_mod2[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                
                                mean_var_mod3[f'{julian_d}']=[]
                                mean_var_mod3[f'{julian_d}'].append({f'{_date}':bn.nanmean(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})
                            
                        except IndexError:
                            pass                    

                        #Index error exists because there is no data at the end of the file
                        #That meets the current idx restrictions
        
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
                                
                                if idx == 0:
                                    mean_var_mod0[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod1[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod2[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                    mean_var_mod3[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=idx).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   

                                
                                elif idx % 7 == 0:
                                    try:
                                        mean_var_mod0[f'{julian_d}'].append({f'{file[day_s:day_e]}':np.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod1[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod2[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                        mean_var_mod3[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   

                                    #Some shouldn't/can't be appended to dictionary because they are useless
                                    except KeyError:
                                        pass
                                    except IndexError:
                                        pass

                                    
                    def find_anomaly(mean_var_modN,mod_number,anomaly_spread):
                        out_dict = {}
                        out_mean = {}
                        '''This will return all the data for each file for averages
                        because of day of year, we have to subset the +/- 42 days in 
                        a different way'''

                        
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
                            #Week 3 & 4 lead
                        
                        out_mean_final = out_mean
                        '''Find new means for weekly averages first
                        Averages are from all of the individual observations'''
                        out_mean_final['3.4'] = ((out_mean['3']+ out_mean['4'])/2).astype('float32')
                        out_mean_final['3.5'] = ((out_mean['3']+ out_mean['4'] + out_mean_final['5']) /3).astype('float32')
                        out_mean_final['3.6'] = ((out_mean['3']+ out_mean['4'] + out_mean['5']+ out_mean['6']) /4).astype('float32')
                        out_mean_final['4.6'] = ((out_mean['4']+ out_mean['5'] + out_mean['6']) /3).astype('float32')
                        
                        #get week lead first from Summation_ETo_modN, then subtract the new mean anomaly
                        
                        def create_multi_weekly_anom (out_dict):
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
                               
                                
                        out_dict = create_multi_weekly_anom (out_dict=out_dict)
                        '''For some reason, the last weekly average lead doesn't work, just type manually'''
                        anom_list = []
                        # elif week_num == '4&5&6':
                        for i in range(len(mean_var_modN[f'{julian_d}'])):
                           anom_list.append(((list(mean_var_modN[list(mean_var_modN.keys())[4-1]][i].values())[0] + \
                                list(mean_var_modN[list(mean_var_modN.keys())[5-1]][i].values())[0] + \
                                             list(mean_var_modN[list(mean_var_modN.keys())[6-1]][i].values())[0])/3) - out_mean_final['4.6'])
                           
                           out_dict['4.6'] = anom_list 
                           

                        return (out_dict,out_mean_final)
                    
                    anom_mod0,mean_mod0 = find_anomaly(mean_var_mod0,mod_number=0,anomaly_spread=anomaly_spread)
                    anom_mod1,mean_mod1 = find_anomaly(mean_var_mod1,mod_number=1,anomaly_spread=anomaly_spread)
                    anom_mod2,mean_mod2 = find_anomaly(mean_var_mod2,mod_number=2,anomaly_spread=anomaly_spread)
                    anom_mod3,mean_mod3 = find_anomaly(mean_var_mod3,mod_number=3,anomaly_spread=anomaly_spread)
                    
                    add_dates = [i[day_s:] for i in dates_to_keep]
                    anom_mod0['init_dates'] = add_dates
                    anom_mod1['init_dates'] = add_dates
                    anom_mod2['init_dates'] = add_dates
                    anom_mod3['init_dates'] = add_dates

                    return(anom_mod0,anom_mod1,anom_mod2,anom_mod3,dates_to_keep, \
                           mean_mod0,mean_mod1,mean_mod2,mean_mod3)
                
                #Contains the julian day value of the current file ETo_{_date} and the 7-day summation
                anom_mod0,anom_mod1,anom_mod2,anom_mod3,file_dates_to_keep, \
                    mean_mod0,mean_mod1,mean_mod2,mean_mod3= dict1_subx2()
                
                #Create weekly leads all in one file
                '''We have created weekly leads for several weeks (see below)
                out_lead = np.array([0,1,2,3,4,5,6,3.4,3.5,3.6,4.6])
                
                Weeks: 1,2,3,4. Then we take averages of multiple weeks including:
                    weeks 3&4, 3&4&5, 3&4&5&6, 4&5&6.
                    
                    Total of 8 different leads
                '''
                
                #Add all date files to it's own dictionary to quickly add to new files
                #This code will add each lead time to it's correct file (verified by my own eyes)
                def improve_anomaly_dictionary( anom_modN,  file_dates_to_keep):
                    #first add a dictinary with the start date as the key
                    mod_out = {}
                    for idx in range(len(file_dates_to_keep)):
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
                
                out_mod0 = improve_anomaly_dictionary(anom_mod0, file_dates_to_keep)
                out_mod1 = improve_anomaly_dictionary(anom_mod1, file_dates_to_keep)
                out_mod2 = improve_anomaly_dictionary(anom_mod2, file_dates_to_keep)
                out_mod3 = improve_anomaly_dictionary(anom_mod3, file_dates_to_keep)


                '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                we can append to aall EDDI files'''
                
                #It would be best to first open 1 file and append all possible values from that one file. Then move onto the next file.
            
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
                        
                        #Now add to file by lead week (or average)
                        for anom_lead in range(len(out_mod0[init_file_day])):

                            #Add data to netcdf file
                            file_open[var_n][0,0,anom_lead,i_Y,i_X] = out_mod0[init_file_day][anom_lead]
                            file_open[var_n][0,1,anom_lead,i_Y,i_X] = out_mod1[init_file_day][anom_lead]
                            file_open[var_n][0,2,anom_lead,i_Y,i_X] = out_mod2[init_file_day][anom_lead]
                            file_open[var_n][0,3,anom_lead,i_Y,i_X] = out_mod3[init_file_day][anom_lead]
                            
                            
                            file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                            file_open.close()

                add_anomaly_to_nc_file(out_mod0, out_mod1,out_mod2,out_mod3)
                
                def add_mean_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep):

                    
                    for idx_,lead in enumerate(mean_mod0):
                        #We have the mean value for each year/lead time/model
                        #Add to only 1 file because its the mean of all years
                                                
                        if var == 'RZSM':
                            day_s = -14
                            fileOut = "{}/{}_mean_{}".format(new_dir_mean,var,file_dates_to_keep[0][day_s:])

                        elif var == 'ETo':
                            day_s = -13
                            fileOut = "{}/{}_mean_{}4".format(new_dir_mean,var,file_dates_to_keep[0][day_s:])

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

                        
                add_mean_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep)
                
                
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
        make_subX_anomaly(_date=_date,var=var, HP_conus_mask = HP_conus_mask, anomaly_spread=42)



# count=0
# for _date in init_date_list[start_:end_]:    
#     if start_ == 50 and _date == '2000-01-10':
#         break
#     else:
#         RZSM_anomaly(start_,end_,init_date_list, _date,var)


