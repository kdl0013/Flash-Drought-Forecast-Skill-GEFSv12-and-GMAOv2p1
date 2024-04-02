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

#%%
# # Test for 1 step size and model
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
start_ = int('0')
end_ = start_ + int('40')
model_NAM1 = 'EMC'
var = 'RZSM'


if var == "RZSM":
    varname='soilw_bgrnd'


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
    _date=init_date_list[0]
    i_X,i_Y=10,10
    anomaly_spread=42
    print(f'Calculating {var} anomaly on SubX for {_date} and saving as .nc4 in {anom_dir}.') 

    #Open all files for faster processing
    if var == 'RZSM':
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}_SubX*.nc4', concat_dim=['S'], combine='nested')
        #Use subx_all to create the full distribution of all data
        # subx_all = xr.open_mfdataset(f'{home_dir}/{varname}*.nc4', combine='by_coords').persist()
        # subx_all = xr.open_mfdataset(f'{home_dir}/{varname}*.nc4', concat_dim=['S'], combine='nested').persist() #works, but lead isn't right
        # subx_all = xr.open_mfdataset(f'{home_dir}/{varname}*.nc4', concat_dim=['S'], combine='nested', chunks={"S": 1},parallel='True').persist()
        subx_all = xr.open_mfdataset(f'{home_dir}/{varname}*.nc4', concat_dim=['S'], combine='nested',parallel='True') #load first

        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{varname}_EMC_{_date}.nc4', engine='netcdf4')

    elif var == 'ETo':
        #Open all files for faster processing (for EDDI and ETo anomaly)
        # subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested')
        subx_all = xr.open_mfdataset(f'{home_dir}/{var}*.nc', concat_dim=['S'], combine='nested').persist()

        #Process indivdual files for output
        subx_out = xr.open_dataset(f'{home_dir}/{var}_{_date}.nc')
        var_name = var
        
    #Restrict lead julian date domain, load into memory early
    julian_d = pd.to_datetime(_date).dayofyear
    
    
    if int(julian_d) <= anomaly_spread:
        subtract_ = int(julian_d)-anomaly_spread
        
        sub_small1 = subx_all[f'{var}'][::,:,:,:].sel(lead=slice(366+subtract_,366))
        sub_small2 = subx_all[f'{var}'][:,:,:,:,:].sel(lead=slice(1,int(julian_d+anomaly_spread)+subx_out.lead.shape[0]))#Add extra lead days
        subx_all = xr.concat([sub_small1,sub_small2],dim='lead').to_dataset().persist()

    elif julian_d >= (366-anomaly_spread):
        add_ = 366-julian_d
        diff_ = anomaly_spread - add_
        sub_small1 = subx_all[f'{var}'][:,:,:,:,:].sel(lead=slice(julian_d,julian_d+add_)) #get dates before new julian_day
        sub_small2 = subx_all[f'{var}'][:,:,:,:,:].sel(lead=slice(julian_d-anomaly_spread,julian_d))
        sub_small3 = subx_all[f'{var}'][:,:,:,:,:].sel(lead=slice(1,diff_+subx_out.lead.shape[0]))
        subx_all = xr.concat([sub_small1,sub_small2,sub_small3],dim='lead').to_dataset().persist()
    else:
        subx_all = subx_all[f'{var}'][:,:,:,:,:].sel(lead=slice(julian_d-anomaly_spread,julian_d+anomaly_spread+subx_out.lead.shape[0])).to_dataset().persist()

    for i_Y in range(subx_out[f'{var}'].shape[3]):
        for i_X in range(subx_out[f'{var}'].shape[4]):
            if _date == init_date_list[0]:
                print(f'Working on lat {i_Y} and lon {i_X}')

            '''Initially, I just wanted to use the SMERGE files as a mask, but it
            does unneccesary computations for a lot of grid cells
            
            #(np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) !=1) or ((i_X == 38 and i_Y == 6) or (i_X ==38 and i_Y ==7))
            #only work on grid cells with values like SMERGE. But just use the USDM Conus mask. Values 1-6 indicate a region with USDM mask.
            
            '''
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

                def find_anomaly_with_weekly_mean(all_mean_var_mod,anomaly_spread=42):
                    out_dict = {}
                    out_mean = {}
                    '''This will return all the data for each file for averages
                    because of day of year, we have to subset the +/- 42 days in 
                    a different way due to how the files day of year is setup'''
                   
                    for mod in all_mean_var_mod.keys():
                        
                        out_mean[mod] = []
                        out_dict[mod] = []
                        for k, julian_d in enumerate((all_mean_var_mod[mod])):
                            julian_d = int(list(julian_d)[0])
                            '''Grab all days with 42 days of day of year, need
                            this weird approach because of slicing issues'''
                            if int(julian_d) <= anomaly_spread:
                                subtract_ = int(julian_d)-anomaly_spread
                                sub_small1 = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(366+subtract_,366))
                                sub_small2 = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(1,int(julian_d+anomaly_spread)))
                                sub_out = xr.concat([sub_small1,sub_small2],dim='lead')

                            elif julian_d >= (366-anomaly_spread):
                                add_ = 366-julian_d
                                diff_ = anomaly_spread - add_
                                sub_small1 = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(julian_d,julian_d+add_)) #get dates before new julian_day
                                sub_small2 = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(julian_d-anomaly_spread,julian_d))
                                sub_small3 = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(1,diff_))
                                sub_out = xr.concat([sub_small1,sub_small2,sub_small3],dim='lead')
                            else:
                                sub_out = subx_all[f'{var}'][:,int(mod),:,i_Y,i_X].sel(lead=slice(julian_d-anomaly_spread,julian_d+anomaly_spread)) 
                
                            '''for some reason, model 3 with ETo has some inf values
                            that are messing up anomaly calculation'''
                            
                            #TODO: This is the bottleneck of the program
                            sub_out_vals = sub_out.values
                            #######################################################
                            sub_out_vals = sub_out_vals.flatten()
                            sub_out_vals[sub_out_vals > 1e06]  = np.nan
                            
                            #now append mean value to dictionary with 
                            mean_value = np.nanmean(sub_out_vals)
                            out_mean[mod].append({str(julian_d):mean_value})
                            # julian_d_str = str(julian_d)
                            # out_dict[mod].append({list(julian_d.keys())[0]:list(julian_d.values())[0] - list(out_mean[mod][0].values())[0]})
                                
                        #subtract mean to make anomaly
                        # julian_d_str = str(julian_d)
                        for idx,julian in enumerate(all_mean_var_mod[mod]):
                            # print(idx,julian)
                            out_dict[mod].append({list(julian.keys())[0]:list(julian.values())[0] - list(out_mean[mod][idx].values())[0]})
                            
                                # 
                                
                                # try:
                                #     value_ = list(mean_var_modN[f'{julian_d}'][i].values())[0]
                                #     anomaly_list.append(value_ - mean_value)
                                # except AttributeError:
                                #     value_=mean_var_modN[f'{julian_d}'][i]
                                #     anomaly_list.append(value_ - mean_value)
                    return(out_mean,out_dict)
                
                model_mean_vals, model_anomaly_values =find_anomaly_with_weekly_mean(all_mean_var_mod,anomaly_spread=42)

                         
                def add_mean_to_nc_file(model_mean_vals):
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

