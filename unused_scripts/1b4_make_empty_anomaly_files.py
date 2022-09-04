o#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove unneeded S dimension in SubX files.

Make empty .nc anomaly files to store weekly leads of:
    All lead values for each model.
    Then weekly averages (may increase skill by reducing noise in data)
    3&4-wk
    3&4&5-wk
    3&4&5&6-wk
    4&5&6-wk

@author: kdl
"""

import xarray as xr
import numpy as np
import os
from glob import glob
import pandas as pd
import datetime as dt
from datetime import timedelta


dir1 = 'main_dir'
mod = 'model_name'
var = 'variables'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod='EMC'
# var='RZSM'

home_dir = f'{dir1}/Data/SubX/{mod}'
out_dir = f'{dir1}/Data/SubX/{mod}/anomaly'
new_acc_dir = f'{out_dir}/mean_for_ACC'
MME_dir = f'{out_dir}/MME'
mme_mean_dir = f'{MME_dir}/mean_MME'

os.system(f'mkdir -p {new_acc_dir}')
os.system(f'mkdir -p {mme_mean_dir}')

script_dir = f'{dir1}/Scripts'
os.chdir(home_dir)

def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'tas')

#For anomaly calculation 
extra_lead_average = [3.4,3.5,3.6,4.6]

if mod=='GMAO' or mod=='RSMAS':
    out_lead = list(np.arange(45)) #45 day lead forecast
    #add extra days for the average
    for idx,val in enumerate(extra_lead_average):
        out_lead.append(extra_lead_average[idx])
elif mod=='EMC':
    out_lead = list(np.arange(35)) #45 day lead forecast
    #add extra days for the average
    for idx,val in enumerate(extra_lead_average):
        out_lead.append(extra_lead_average[idx])
elif mod=='ESRL':
    out_lead = list(np.arange(32)) #45 day lead forecast
    #add extra days for the average
    for idx,val in enumerate(extra_lead_average):
        out_lead.append(extra_lead_average[idx])


#%%'''Make empty files to keep track of anomaly calculations'''

#Make a completed list file for EDDI, add new names to next code block
new_eddi = f'EDDI_completed_nc_{mod}.txt'

new_eto_anom = f'ETo_completed_mean_for_anomaly_{mod}.txt'
# new_eto_mean= f'ETo_completed_mean_nc_{mod}.txt'

new_rzsm = f'RZSM_completed_mean_{mod}_for_anomaly_.txt'
# new_rzsm_mean = f'RZSM_completed_mean_nc_{mod}.txt'

# #Make a completed list file for EDDI, add new names to next code block
# new_eddi1 = f'EDDI_MME_completed_nc_{mod}.txt'

# new_eto_anom1 = f'ETo_MME_completed_anomaly_nc_{mod}.txt'
# # new_eto_mean= f'ETo_completed_mean_nc_{mod}.txt'

# new_rzsm1 = f'RZSM_MME_completed_anomaly_nc_{mod}.txt'

# new_rzsm_mean = f'RZSM_completed_mean_nc_{mod}.txt'
# [new_eddi,new_eto_anom,new_rzsm,new_eddi1,new_eto_anom1,new_rzsm1]
#Create a new file for each index to keep track of what has been completed 
#since this is a pretty slow process
for name in [new_eddi,new_eto_anom,new_rzsm]:
    try:
        completed_dates = np.loadtxt(f'{script_dir}/{name}',dtype='str')         
    except OSError:
        os.system(f'touch {script_dir}/{name}')
        name_index = f'{name}'.split('_')[0]
        os.system(f'echo "Completed {name_index}" > {script_dir}/{name}')

#%%
T_FILE = xr.open_dataset('ETo_2000-01-15.nc')

def make_empty_mean_nc_files(init_date_list,T_FILE, var,out_lead):
    print(f'Making empty .nc files for {var} and saving into {new_acc_dir} [if applicable].')
    # init_date_list[0][5:]
    #Return date list for a single year with all possible day values from other files
    tmp = [i[5:] for i in init_date_list]
    init_date_out = []

    [init_date_out.append(date) for date in tmp if date not in init_date_out]    
    
    for _date in init_date_out:
        #Make the file the year 2000
        model_dirs = {}
        # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
        if var == 'EDDI':
            filename = 'EDDI_mean'
            desc = 'Evaporative Demand Drought Index mean'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                
        elif var == 'RZSM':
            filename = f'{var}_mean'
            desc = f'RZSM mean SubX {mod} model.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
        elif var == 'ETo':
            filename = f'{var}_mean'
            desc = f'ETo mean SubX {mod} model.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
           
            
        #don't remake empty files
        try:
            xr.open_dataset(f'{new_acc_dir}/{filename}_2000-{_date}.nc4')
            
        except FileNotFoundError:
            template = xr.open_dataset(f'ETo_{_date}.nc')
            #initialize empty file
            template = xr.zeros_like(template)
            template_out = template.ETo[:,:,0:len(out_lead),:,:]
        
            
            if var == 'RZSM':
                var_final = xr.Dataset(
                    data_vars = dict(
                        RZSM_mean = (['S','model','lead','Y','X'], template_out.data),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                
            elif var == 'EDDI':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        EDDI_mean = (['S','model','lead','Y','X'], template_out.data),
                    ),   
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                

                
            elif var == 'ETo':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        ETo_mean = (['S','model','lead','Y','X'], template_out.data),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values,
                        model = template.model.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                

            var_final.to_netcdf(path = f'{new_acc_dir}/{filename}_2000-{_date}.nc4')
            #compress so that when I re-write the file, it is quicker
            # '''But this doesn't work after the file is re-read and re-saved'''
            # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
            # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
            var_final.close()
            
make_empty_mean_nc_files(init_date_list,T_FILE, var,out_lead)

#%%
# #I brought this function out of loop, no need to repeat it more than once
# if mod == 'GMAO':
#     T_FILE = xr.open_dataset('ETo_2000-01-15.nc')

# def make_empty_nc_files(init_date_list,T_FILE, var,out_lead):
#     print(f'Making empty .nc files for {var}.')
#     # count=0
#     # _date=init_date_list[0]

#     for _date in init_date_list:
       
#         model_dirs = {}
#         # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
#         if var == 'EDDI':
#             filename = 'EDDI'
#             desc = 'Evaporative Demand Drought Index. Calculated by lead week (1-6) over all 15 years of dataset. Weekly leads of \
#                 1,2,3,4, mean 3&4, mean 3&4&5, mean 4&5&6'''
#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                
#         elif var == 'RZSM':
#             filename = f'{var}_anomaly'
#             desc = f'RZSM anomaly SubX {mod} model. Calculated by lead week (1-6) over all 15 years of dataset. Weekly leads of \
#                 1,2,3,4, mean 3&4, mean 3&4&5, mean 4&5&6'''
               
#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
#         elif var == 'ETo':
#             filename = f'{var}_anomaly'
#             desc = f'ETo anomaly SubX {mod} model. Calculated by lead week (1-6) over all 15 years of dataset. Weekly leads of \
#                 1,2,3,4, mean 3&4, mean 3&4&5, mean 4&5&6'''

#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
        
#         # try:
#         #     np.load(f'{home_dir}/julian_lead_{_date}.npy')
#         # except FileNotFoundError:
#         #     np.save(f'{home_dir}/julian_lead_{_date}.npy',day_julian_b) 
            
            
#         #don't remake empty files
        
#         try:
#             xr.open_dataset(f'{out_dir}/{filename}_{_date}.nc4')

#         except FileNotFoundError:
#             template = xr.open_dataset(f'ETo_{_date}.nc')
#             #initialize empty file
#             template = xr.zeros_like(template)
            
#             template_out = template.ETo[:,:,0:len(out_lead),:,:]
            
#             if var == 'RZSM':
#                 var_final = xr.Dataset(
#                     data_vars = dict(
#                         RZSM_anom = (['S','model','lead','Y','X'], template_out.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values,
#                         model = template.model.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc}'),
#                 )             
                

#             elif var == 'EDDI':
#                 #save file as EDDI as netcdf
#                 var_final = xr.Dataset(
#                     data_vars = dict(
#                         EDDI = (['S','model','lead','Y','X'], template_out.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values,
#                         model = template.model.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc}'),
#                 )                

                
#             elif var == 'ETo':
#                 #save file as EDDI as netcdf
#                 var_final = xr.Dataset(
#                     data_vars = dict(
#                         ETo_anom = (['S','model','lead','Y','X'], template_out.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values,
#                         model = template.model.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc}'),
#                 )                    
    
#             var_final.to_netcdf(path = f'{out_dir}/{filename}_{_date}.nc4')
#             #compress so that when I re-write the file, it is quicker
#             # '''But this doesn't work after the file is re-read and re-saved'''
#             # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
#             # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
#             var_final.close()

   
# #Create EDDI files           
# make_empty_nc_files(init_date_list = init_date_list,T_FILE = T_FILE, var=var,out_lead=out_lead)

#%%
# def make_empty_multi_member_ensemble_files(init_date_list = init_date_list,T_FILE = T_FILE, var=var,out_lead=out_lead):
    
#     for _date in init_date_list:
       
#         model_dirs = {}
#         # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
#         if var == 'EDDI':
#             filename = 'EDDI_MME'
#             desc = 'Multi-member ensemble mean. Evaporative Demand Drought Index.'''
#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
#         elif var == 'RZSM':
#             filename = f'{var}_MME_anomaly'
#             desc2 = f'Multi-member ensemble mean. RZSM SubX {mod} model.'''                
#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
#         elif var == 'ETo':
#             filename = f'{var}_MME_anomaly'
#             desc2 = f'Multi-member ensemble mean. ETo anomaly SubX {mod} model.'''
#             for model in [0,1,2,3]:
#                 model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
           
#         #don't remake empty files
#         try:
#             xr.open_dataset(f'{MME_dir}/{filename}_{_date}.nc4')
#         except FileNotFoundError:
#             template = xr.open_dataset(f'ETo_{_date}.nc')
#             #initialize empty file
#             template = xr.zeros_like(template)
            
#             template_out = template.ETo[:,:,0:len(out_lead),:,:]
#             template_out_MME = template.ETo[:,0,0:len(out_lead),:,:]
            
#             if var == 'RZSM':
#                 var_final_MME = xr.Dataset(
#                     data_vars = dict(
#                         RZSM_MME_anom = (['S','lead','Y','X'], template_out_MME.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc2}'),
#                 )             
#             elif var == 'EDDI':
#                 #save file as EDDI as netcdf
#                 var_final_MME = xr.Dataset(
#                     data_vars = dict(
#                         EDDI_MME = (['S','lead','Y','X'], template_out_MME.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc}'),
#                 )                

                
#             elif var == 'ETo':
#                 #save file as ETo as netcdf
#                 var_final_MME = xr.Dataset(
#                     data_vars = dict(
#                         ETo_MME_anom = (['S','lead','Y','X'], template_out_MME.data),
#                     ),
#                     coords = dict(
#                         X = template.X.values,
#                         Y = template.Y.values,
#                         lead = out_lead,
#                         S = template.S.values
#                     ),
#                     attrs = dict(
#                         Description = f'{desc2}'),
#                 )             
#             var_final_MME.to_netcdf(path = f'{MME_dir}/{filename}_{_date}.nc4')
#             #compress so that when I re-write the file, it is quicker
#             # '''But this doesn't work after the file is re-read and re-saved'''
#             # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
#             # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
#             var_final_MME.close()

# make_empty_multi_member_ensemble_files(init_date_list = init_date_list,T_FILE = T_FILE, var=var,out_lead=out_lead)
#%% Now make empty mean files for ACC
#%%

#%%

def make_empty_MME_mean_nc_files(init_date_list,T_FILE, var,out_lead):
    print(f'Making empty .nc files for {var} and saving into {mme_mean_dir} [if applicable].')


    for _date in init_date_list[0:73]:
       
        model_dirs = {}
        # S_values = [pd.to_datetime(zz)+ dt.timedelta(days=1),pd.to_datetime(zz)]
        if var == 'EDDI':
            filename = 'EDDI_MME_mean'
            desc = 'Evaporative Demand Drought Index multi member ensemble mean.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/EDDI_mod{model}'
                
        elif var == 'RZSM':
            filename = f'{var}_MME_mean'
            desc = f'RZSM multi member ensemble mean SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
            
        elif var == 'ETo':
            filename = f'{var}_MME_mean'
            desc = f'ETo multi member ensemble mean SubX {mod} model. Calculated by lead week (1-7) over all 15 years of dataset.'
            for model in [0,1,2,3]:
                model_dirs[f'Model {model}'] = f'{home_dir}/{var}_anomaly_mod{model}'
           
            
        #don't remake empty files
        try:
            xr.open_dataset(f'{new_acc_dir}/{filename}_{_date}.nc4')
            
        except FileNotFoundError:
            template = xr.open_dataset(f'ETo_{_date}.nc')
            #initialize empty file
            template = xr.zeros_like(template)
            template_out = template.ETo[:,:,0:len(out_lead),:,:]
            template_out = template.ETo[:,0,0:len(out_lead),:,:] #combine all models
            
            if var == 'RZSM':
                var_final = xr.Dataset(
                    data_vars = dict(
                        RZSM_MME_mean = (['S','lead','Y','X'], template_out.data),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                
            elif var == 'EDDI':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        EDDI_MME_mean = (['S','lead','Y','X'], template_out.data),
                    ),   
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                

                
            elif var == 'ETo':
                #save file as EDDI as netcdf
                var_final = xr.Dataset(
                    data_vars = dict(
                        ETo_MME_mean = (['S','lead','Y','X'], template_out.data),
                    ),
                    coords = dict(
                        X = template.X.values,
                        Y = template.Y.values,
                        lead = out_lead,
                        S = template.S.values
                    ),
                    attrs = dict(
                        Description = f'{desc}'),
                )                    
                

            var_final.to_netcdf(path = f'{mme_mean_dir}/{filename}_{_date}.nc4')
            #compress so that when I re-write the file, it is quicker
            # '''But this doesn't work after the file is re-read and re-saved'''
            # os.system(f'ncks -4 -L 1 {home_dir}/{filename}_{_date}.nc4 {home_dir}/{filename}_{_date}_test.nc4')
            # os.system(f'mv {home_dir}/{filename}_{_date}_test.nc4 {home_dir}/{filename}_{_date}.nc4')
            
            var_final.close()
            
make_empty_MME_mean_nc_files(init_date_list,T_FILE, var,out_lead)
