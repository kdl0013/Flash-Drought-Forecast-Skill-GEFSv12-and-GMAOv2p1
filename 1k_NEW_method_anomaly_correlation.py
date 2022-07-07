#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Steps:
    
    1.) re use 1c_make_anomaly code to re-get the actual forecasts and then just subtract
    the climatology mean (already created (but slow process initially)) from each date.
    This will allow me to create the variance which is the bottom portion of the 
    ACC formula.
    2.) Then open up the observations (with the associated observations mean file)
    and 

Source:
https://metclim.ucd.ie/wp-content/uploads/2017/07/DeterministicSkillScore.pdf
def ACC(FC,OBS,CL):
    top = np.mean((FC-CL)*(OBS-CL))
    bottom = np.sqrt(np.mean((FC-CL)**2)*np.mean((OBS-CL)**2))
    ACC = top/bottom
    return ACC

@author: kdl
"""


import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.style.use('seaborn-talk')

import xarray as xr
import numpy as np
from numba import njit, prange
import numpy.ma as ma
import os
import datetime as dt
import pandas as pd
from glob import glob
from scipy.stats import rankdata
import xskillscore as xs
import sys
import dask
from climpred import HindcastEnsemble
import climpred
import cartopy.crs as ccrs
import cartopy as cart
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl


# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'

os.chdir(f'{dir1}/Data/SubX/{model_NAM1}')

#TODO change later
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod = 'GMAO'

dir1 = 'main_dir'
model_NAM1 = 'model_name'
subX_dir = f'{dir1}/Data/SubX/{model_NAM1}/' #where subX model data lies 

subdir_for_mean = f'{subX_dir}/mean_for_ACC'
mask_path = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4'

#observations
obs_eto_path = f'{dir1}/Data/gridMET/ETo_SubX_values/'
obs_rzsm_path = f'{dir1}/Data/SMERGE_SM/SM_SubX_values/'

obs_eto_mean = xr.open_dataset(f'{dir1}/Data/gridMET/ETo_anomaly_gridMET_mean.nc')
obs_rzsm_mean = xr.open_dataset(f'{dir1}/Data/SMERGE_SM/Raw_data/RZSM_anomaly_SMERGE_mean.nc')

output_nc_dir = f'{subX_dir}/skill_assessments/' #for skill assessment .nc files
output_image_dir = f'{dir1}/Outputs/pearson_correlation/'
output_season_dir = f'{dir1}/Outputs/pearson_correlation/seasonal_skill/'


var_list='mrso'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{subX_dir}/{var_list}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list()    

#%%

# _date=init_date_list[0]
print(f'Calculating {var} ACC on SubX for {_date} and saving as .nc4 in {home_dir}.') 
#Used for eliminating iterating over grid cells that don't matter
smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
smerge_file_julian = smerge_file.copy()
  
#Open all files for faster processing
if var == 'RZSM':
    #Process indivdual files for output
    subx_out = xr.open_dataset(f'{home_dir}/{var}_anomaly_{_date}.nc4')
    subx_out.close()
    var_name='RZSM_anom'


elif var == 'ETo':
    subx_out = xr.open_dataset(f'{home_dir}/{var}_anomaly_{_date}.nc4')
    subx_out.close()
    var_name='ETo_anom'

#Convert to julian day for processing
smerge_day = pd.to_datetime(smerge_file.CCI_ano.time.values)
smerge_julian = [i.timetuple().tm_yday for i in smerge_day]

smerge_file_julian= smerge_file_julian.assign_coords({'time': smerge_julian})

for i_Y in range(subx_out[f'{var_name}'].shape[3]):
    for i_X in range(subx_out[f'{var_name}'].shape[4]):
        if _date == '1999-01-10':
            print(f'Working on lat {i_Y} and lon {i_X}')

        #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
        #only work on grid cells with values like SMERGE
        
        if var == 'RZSM':
            true_or_false = (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1)
        elif var == 'ETo':
            true_or_false = (np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
        
        #This is for the mask of CONUS, don't do extra unneeded calculations
        if true_or_false:
            
            def dict1_subx2():
                
                #append 7-day summation from all files to a new dictionary
                #Technically, this is RZSM, but its easier to keep as ETo for EDDI, RZSM, and ETo
                summation_ETo_mod0 = {}
                summation_ETo_mod1 = {}
                summation_ETo_mod2 = {}
                summation_ETo_mod3 = {}
                
                for idx,julian_d in enumerate(subx_out.lead.values):
                    #You must julian_d + week_lead because with RZSM you need 7-day looking backwards into the past. Must have 7 values.
                    #Choose just model = 0 because we just need to know if there are 7-days total in any model
                    # print(idx)

                    try:
                        if idx % 7 == 0 and idx !=0:
                            summation_ETo_mod0[f'{julian_d}']=[]
                            summation_ETo_mod0[f'{julian_d}'].append({f'{_date}':bn.nansum(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                            
                            summation_ETo_mod1[f'{julian_d}']=[]
                            summation_ETo_mod1[f'{julian_d}'].append({f'{_date}':bn.nansum(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                            
                            summation_ETo_mod2[f'{julian_d}']=[]
                            summation_ETo_mod2[f'{julian_d}'].append({f'{_date}':bn.nansum(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                            
                            summation_ETo_mod3[f'{julian_d}']=[]
                            summation_ETo_mod3[f'{julian_d}'].append({f'{_date}':bn.nansum(subx_out[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})
                        
                    except IndexError:
                        pass                    

                    #Index error exists because there is no data at the end of the file
                    #That meets the current idx restrictions
    
                '''7-day mean of RZSM by lead date for the first opened file:
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
                            if idx % 7 == 0 and idx != 0:
                                try:
                                    summation_ETo_mod0[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=0, X=i_X, Y=i_Y).values)})   
                                    summation_ETo_mod1[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=1, X=i_X, Y=i_Y).values)})   
                                    summation_ETo_mod2[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=2, X=i_X, Y=i_Y).values)})   
                                    summation_ETo_mod3[f'{julian_d}'].append({f'{file[day_s:day_e]}':bn.nanmean(open_f[f'{var_name}'].isel(lead=slice(idx-7,idx)).isel(S=0, model=3, X=i_X, Y=i_Y).values)})   

                                #Some shouldn't/can't be appended to dictionary because they are useless
                                except KeyError:
                                    pass
                                except IndexError:
                                    pass

                                
                def find_anomaly(summation_ETo_modN,mod_number,anomaly_spread):
                    out_dict = {}
                    out_mean = {}
                    '''This will return all the data for each file for averages'''
                    for k, julian_d in enumerate(summation_ETo_modN):
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
 
                        mean_value = bn.nanmean(sub_out)
                        out_mean[f'{julian_d}'] = mean_value
                        #subtract mean
                        anomaly_list = []
                        julian_d = str(julian_d)
                        for i in range(len(summation_ETo_modN[f'{julian_d}'])):
                            try:
                                value_ = list(summation_ETo_modN[f'{julian_d}'][i].values())[0]
                                anomaly_list.append(value_ - mean_value)
                            except AttributeError:
                                value_=summation_ETo_modN[f'{julian_d}'][i]
                                anomaly_list.append(value_ - mean_value)


                        out_dict[f'{julian_d}']=anomaly_list
                        
                    return (out_dict,out_mean)
                
                anom_mod0,mean_mod0 = find_anomaly(summation_ETo_mod0,mod_number=0,anomaly_spread=anomaly_spread)
                anom_mod1,mean_mod1 = find_anomaly(summation_ETo_mod1,mod_number=1,anomaly_spread=anomaly_spread)
                anom_mod2,mean_mod2 = find_anomaly(summation_ETo_mod2,mod_number=2,anomaly_spread=anomaly_spread)
                anom_mod3,mean_mod3 = find_anomaly(summation_ETo_mod3,mod_number=3,anomaly_spread=anomaly_spread)
                
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
            def improve_anomaly_dictionary( anom_modN,  file_dates_to_keep):
                #first add a dictinary with the start date as the key
                mod_out = {}
                for idx in range(len(file_dates_to_keep)):
                    mod_out[anom_modN['init_dates'][idx]] = []
                    
                    all_values = []
                    for init_day,values in anom_modN.items():

                        if init_day=='init_dates':
                            break
                        else:
                            all_values.append(anom_modN[init_day][idx])
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
                        fileOut = ("{}_anomaly_{}4".format(var,init_file_day))
                    elif var == 'RZSM':
                        fileOut = ("{}_anomaly_{}".format(var,init_file_day))
                        
                    var_n = f'{var}_anom'
                    file_open = xr.open_dataset(fileOut)
                    file_open.close()
                    
                    #Now add to file by lead week
                    lead_week=1
                    for anom_vals in range(len(out_mod0[init_file_day])):

                        #Add data to netcdf file
                        file_open[var_n][0,0,lead_week*7,i_Y,i_X] = out_mod0[init_file_day][anom_vals]
                        file_open[var_n][0,1,lead_week*7,i_Y,i_X] = out_mod1[init_file_day][anom_vals]
                        file_open[var_n][0,2,lead_week*7,i_Y,i_X] = out_mod2[init_file_day][anom_vals]
                        file_open[var_n][0,3,lead_week*7,i_Y,i_X] = out_mod3[init_file_day][anom_vals]
                        
                        
                        file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                        file_open.close()
                        lead_week+=1

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
                    file_open[var_n][0,0,idx_*7,i_Y,i_X] = mean_mod0[f'{lead}']
                    file_open[var_n][0,1,idx_*7,i_Y,i_X] = mean_mod1[f'{lead}']
                    file_open[var_n][0,2,idx_*7,i_Y,i_X] = mean_mod2[f'{lead}']
                    file_open[var_n][0,3,idx_*7,i_Y,i_X] = mean_mod3[f'{lead}']
                    
                    file_open.to_netcdf(path = fileOut, mode ='w', engine='scipy')
                    file_open.close()

                    
            add_mean_to_nc_file(mean_mod0,mean_mod1,mean_mod2,mean_mod3,file_dates_to_keep)
            
            
print(f'Completed date {_date} and saved into {home_dir}.')

#save the dates that were completed to not re-run
os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_anomaly_nc_{model_NAM1}.txt')
os.system(f'echo Completed {_date} >> {script_dir}/{var}_completed_mean_nc_{model_NAM1}.txt')

return()
#%%



# _date=init_date_list[0]
'''Read RZSM_completed_anomaly_nc_.txt file to not have to re-run extra code'''
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



