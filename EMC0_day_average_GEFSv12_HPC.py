#!/usr/bin/env python3

'''Because all GEFSv12 models are in seperate files for each model, I need to save
on space. So combine all models into a single file for each model and for each lead day.

Since there are two files for each variable (days 0-10 and another for days 10-35), 
also need to combine these into a single file for a 35 day forecast.

https://noaa-gefs-retrospective.s3.amazonaws.com/Description_of_reforecast_data.pdf

'''
import os
import datetime as dt
import numpy as np
import xarray as xr
from glob import glob
import pandas as pd
from multiprocessing import Pool


dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'

# #for linux
home_dir = f'{dir1}/Data/SubX/EMC/raw_ensemble_files/'
save_out_dir_GEFS = f'{dir1}/Data/SubX/EMC/'

vars_to_process= [i for i in os.listdir(home_dir)]
# for i in vars_to_process:
#     os.system(f'mkdir -p {out_dir}/{i}')
#%%

def merge_ensemble_members(var):
# for var in vars_to_process:
    # var=vars_to_process[0]
    print(f'Working on variable {var} to merge ensemble members.')

    os.chdir(f'{home_dir}/{var}')
    
    #dates
    #GEFS long-term (multi-ensemble) forecasts are only initialized on Wednesdays
    start_date = dt.date(2000, 1, 1)
    #Actual end date (only data through that period for this website https://noaa-gefs-retrospective.s3.amazonaws.com/index.html#GEFSv12/reforecast/)
    end_date = dt.date(2019, 12, 31)
    dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]
    #from date time, Wednesday is a 2. (Monday is a 0) https://docs.python.org/3/library/datetime.html#datetime.datetime.weekday
    dates = [i for i in dates if i.weekday() ==2]

    '''Need to create a new xarray template with 11 ensemble members, 35 lead days, lat/lon same as GMAO 
    (because file is already regridded to be in the same format as all other observations and GMAO '''


    # def merge_ensemble_members_multiprocess(var):
    #This template was created from the code above
    template_GEFS_initial = np.empty(shape=(1,11,35,27,59))
    lead_splices = ['d10','d35']

    all_possible_ensemble_members = ['c00', 'p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09', 'p10']
    #testing with one at a time
    # vars_to_process = ['soilw_bgrnd']

    #Get the dates of the files
    for _date in dates:
        # _date=dates[0]

        path_out = f"{save_out_dir_GEFS}"
        #The date of the file is when it was intialized
        out_date_create = pd.to_datetime(_date) 
        out_date = f'{out_date_create.year}-{out_date_create.month:02}-{out_date_create.day:02}'
                
        #Make the names the same as the 
        
        if var=='apcp_sfc':
            final_out_name = f'pr_EMC_{out_date}.nc'
        elif var=='tmp_2m':
            final_out_name = f'tas_EMC_{out_date}.nc'
        elif var=='dswrf_sfc':
            final_out_name = f'dswrf_EMC_{out_date}.nc'
        elif var=='dlwrf_sfc':
            final_out_name = f'dlwrf_EMC_{out_date}.nc'
        elif var=='ulwrf_sfc':
            final_out_name = f'ulwrf_EMC_{out_date}.nc'
        elif var=='uswrf_sfc':
            final_out_name = f'uswrf_EMC_{out_date}.nc'
        elif var=='soilw_bgrnd':
            final_out_name = f'RZSM_EMC_{out_date}.nc'
        elif var=='spfh_2m':
            final_out_name = f'huss_EMC_{out_date}.nc'
        elif var=='tmax_2m':
            final_out_name = f'tasmax_EMC_{out_date}.nc'
        elif var=='tmin_2m':
            final_out_name = f'tasmin_EMC_{out_date}.nc'
        elif var=='uflx_sfc':
            final_out_name = f'uas_EMC_{out_date}.nc'
        elif var=='vflx_sfc':
            final_out_name = f'vas_EMC_{out_date}.nc'        
        # final_out_name= f"{var.split('_')[0]}_{out_date}.nc"
        try:
            xr.open_dataset(f"{path_out}/{final_out_name}4")
            #uncomment below to recreate new files
            # t_name='test.nc'
            # xr.open_dataset(f'{t_name}')
        except FileNotFoundError:
            
            template_GEFS_initial[:,:,:,:,:] = np.nan
            # lead_day=0 #keeps up with which index is correct in template_GEFS_initial (resets with each date)

            all_files_d10 = sorted(glob(f'{lead_splices[0]}_{var}_{_date.year}{_date.month:02}{_date.day:02}*.nc4'))
            all_files_d35 = sorted(glob(f'{lead_splices[1]}_{var}_{_date.year}{_date.month:02}{_date.day:02}*.nc4'))

            #some files have doubles (rsync error or from HPC when converting)
            file_len = [len(i) for i in all_files_d10]
            mode = max(set(file_len), key=file_len.count)
            #Replace files
            all_files_d10 = [i for i in all_files_d10 if len(i) == mode]
            all_files_d35 = [i for i in all_files_d35 if len(i) == mode]
               
                
            #TODO:If all realizations are present
            if (len(all_files_d10) == 11) and (len(all_files_d35) == 11):
                #If all possible files are there, then this is the easy code processing to add data to single file
                for ensemble_number,files in enumerate(zip(all_files_d10,all_files_d35)):
             
                    if var != 'soilw_bgrnd':
                        open_d10 = xr.open_dataset(files[0])
                        open_d35 = xr.open_dataset(files[1])
                        var_name = [i for i in list(open_d10.keys()) if 'step' not in i][0]
                    else:
                        open_d10= xr.open_dataset(files[0])
                        open_d35 = xr.open_dataset(files[1])
                        var_name = [i for i in list(open_d10.keys()) if 'step' not in i][0]
                        #TODO: Take the summation of the first 3 soil layers (0-100cm)
                        open_d10 = open_d10[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                        open_d35 = open_d35[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                        

                    #First get the dates of the files
                    '''Take average of first 7 timesteps if d10 file. I have verified
                    this is correct when looking at HPC'''
                    start_ = 0
                    steps = {}
                    for i in range(35):
                        if i ==0:
                            steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values
                            start_+=8 #needed to begin the next index to keep up with proper dates
                        elif i<10:
                            steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values #eight total possible values until last time step
                            start_+= 8 #Need to add one because we don't want to re-index the same day
                        elif i == 10:
                            try:
                                #Need to take from the first file (time 00:00:00), and combine with d35 files
                                s1 = (open_d35[f'{var_name}'][-1,:,:] + open_d35[f'{var_name}'][0,:,:] + \
                                    open_d35[f'{var_name}'][1,:,:] + open_d35[f'{var_name}'][2,:,:]) /4
                                steps[i] = s1
                                start_ = 3 #start count over, 4th file is the new date in d35 files
                            except IndexError:
                                pass
                                #Some ensembles have broken members
                        elif i <=34:
                            steps[i] = open_d35[f'{var_name}'][start_:start_+4,:,:].mean(dim=['step']).values
                            start_+=4
                    #Add to file
                    for step,lead_day in enumerate(steps.keys()):
                        template_GEFS_initial[:,ensemble_number,step,:,:] = steps[lead_day]

            #If all ensembles are missing, do nothing
            elif (len(all_files_d10) == 0 )and (len(all_files_d35) == 0):
                pass
            
            #If there are a differnet number of ensembles between leads
            elif (len(all_files_d10) != 11) or (len(all_files_d35) != 11):

                #Some ensembles are missing, split to get the name of ensemble members
                avail_ensemble_members_d10 = [i.split('_')[-1].split('.')[0] for i in all_files_d10]
                avail_ensemble_members_d35 = [i.split('_')[-1].split('.')[0] for i in all_files_d35]
               
                #if missing only the exact same data
                if len(list(set(avail_ensemble_members_d10).difference(avail_ensemble_members_d35))) == 0:
                    #Find a way to append the missing ensemble files with np.nan
                    
                    for idx,ensemble in enumerate(all_possible_ensemble_members):
                        if ensemble not in avail_ensemble_members_d10:
                            pass
                        else:
                            idx_num = avail_ensemble_members_d10.index(ensemble)
                            open_d10=xr.open_dataset(all_files_d10[idx_num])
                            open_d35=xr.open_dataset(all_files_d35[idx_num])
                            var_name = [i for i in list(open_d10.keys()) if 'step' not in i][0]
                            
                            if var == 'soilw_bgrnd':
                                open_d10 = open_d10[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                                open_d35 = open_d35[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                          
                                #TODO: Take the summation of the first 3 soil layers (0-100cm)
   
                            #First get the dates of the files
                            '''Take average of first 7 timesteps if d10 file. I have verified
                            this is correct when looking at HPC'''
                            start_ = 0
                            steps = {}
                            for i in range(35):
                                if i ==0:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values
                                    start_+=8 #needed to begin the next index to keep up with proper dates
                                elif i<10:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values #eight total possible values until last time step
                                    start_+= 8 #Need to add one because we don't want to re-index the same day
                                elif i == 10:
                                    #Need to take from the first file (time 00:00:00), and combine with d35 files
                                    s1 = (open_d35[f'{var_name}'][-1,:,:] + open_d35[f'{var_name}'][0,:,:] + \
                                        open_d35[f'{var_name}'][1,:,:] + open_d35[f'{var_name}'][2,:,:]) /4
                                    steps[i] = s1
                                    start_ = 3 #start count over, 4th file is the new date in d35 files
                                elif i <=34:
                                    steps[i] = open_d35[f'{var_name}'][start_:start_+4,:,:].mean(dim=['step']).values
                                    start_+=4
                            #Add to file
                            for step,lead_day in enumerate(steps.keys()):
                                template_GEFS_initial[:,idx,step,:,:] = steps[lead_day]
                                
                                
                #if missing different ensemble members -- ONLY 1 MISSING MEMBER
                elif len(list(set(avail_ensemble_members_d10).difference(avail_ensemble_members_d35))) ==1:
                     #Because there are missing ensemble members that are supposed to be aligned,
                     #we must delete those ensemble members
                      
                     missing_member_ = list(set(avail_ensemble_members_d10).difference(avail_ensemble_members_d35))[0]
                     
                     out_10 = []
                     for i in avail_ensemble_members_d10:
                         if missing_member_ not in i:
                             out_10.append(i)
                             
                     out_35 = []
                     for i in avail_ensemble_members_d35:
                         if missing_member_ not in i:
                             out_35.append(i)
                             
                             
                     for idx,ensemble in enumerate(all_possible_ensemble_members):
                         if (ensemble not in out_10) and (ensemble not in out_35):
                             #Value is already in the file, so just pass
                             pass  
                         else:

                            #Current issue. We need to open up the files after removing the missing
                            #ensemble member from one of the lead files
                            
                            all_d10_files = [i for i in all_files_d10 if missing_member_ not in i]
                            all_d35_files = [i for i in all_files_d35 if missing_member_ not in i]
                            #rename briefly to assit
                            d10_pre=[i for i in all_d10_files if ensemble in i][0]
                            d35_pre=[i for i in all_d35_files if ensemble in i][0]
                            
                            open_d10=xr.open_dataset(d10_pre)
                            open_d35=xr.open_dataset(d35_pre)
                            var_name = [i for i in list(open_d10.keys()) if 'step' not in i][0]
                            
                            if var == 'soilw_bgrnd':
                                open_d10 = open_d10[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                                open_d35 = open_d35[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()

                            #First get the dates of the files
                            '''Take average of first 7 timesteps if d10 file. I have verified
                            this is correct when looking at HPC'''
                            start_ = 0
                            steps = {}
                            for i in range(35):
                                if i ==0:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values
                                    start_+=8 #needed to begin the next index to keep up with proper dates
                                elif i<10:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values #eight total possible values until last time step
                                    start_+= 8 #Need to add one because we don't want to re-index the same day
                                elif i == 10:
                                    #Need to take from the first file (time 00:00:00), and combine with d35 files
                                    s1 = (open_d35[f'{var_name}'][-1,:,:] + open_d35[f'{var_name}'][0,:,:] + \
                                        open_d35[f'{var_name}'][1,:,:] + open_d35[f'{var_name}'][2,:,:]) /4
                                    steps[i] = s1
                                    start_ = 3 #start count over, 4th file is the new date in d35 files
                                elif i <=34:
                                    steps[i] = open_d35[f'{var_name}'][start_:start_+4,:,:].mean(dim=['step']).values
                                    start_+=4
                            #Add to file
                            for step,lead_day in enumerate(steps.keys()):
                                template_GEFS_initial[:,idx,step,:,:] = steps[lead_day] 
                                
                #if missing different ensemble members -- MISSING 2 MEMBERS
                elif len(list(set(avail_ensemble_members_d10).difference(avail_ensemble_members_d35))) ==2:
                     #Because there are missing ensemble members that are supposed to be aligned,
                     #we must delete those ensemble members
                      
                     missing_members = list(set(avail_ensemble_members_d10).difference(avail_ensemble_members_d35))
                     
                     #remove missing members
                     out_10 = [i for i in avail_ensemble_members_d10 if i not in missing_members]
                     out_35 = [i for i in avail_ensemble_members_d35 if i not in missing_members]
                             
                     for idx,ensemble in enumerate(all_possible_ensemble_members):
                         if (ensemble not in out_10) and (ensemble not in out_35):
                             #Value is already in the file, so just pass
                             pass  
                         else:
                             
                            #Current issue. We need to open up the files after removing the missing
                            #ensemble member from one of the lead files
                            all_d10_files = [i for i in all_files_d10 if missing_members[0] if missing_members[1] not in i]
                            all_d35_files = [i for i in all_files_d35 if missing_members[0] if missing_members[1] not in i]
                            #rename briefly 
                            d10_pre=[i for i in all_d10_files if ensemble in i][0]
                            d35_pre=[i for i in all_d35_files if ensemble in i][0]
                            
                            open_d10=xr.open_dataset(d10_pre)
                            open_d35=xr.open_dataset(d35_pre)
                            var_name = [i for i in list(open_d10.keys()) if 'step' not in i][0]
                            
                            if var == 'soilw_bgrnd':
                                open_d10 = open_d10[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()
                                open_d35 = open_d35[f'{var_name}'][:,0:2,:,:].sum(dim=['depthBelowLandLayer']).to_dataset()

                            #First get the dates of the files
                            '''Take average of first 7 timesteps if d10 file. I have verified
                            this is correct when looking at HPC'''
                            start_ = 0
                            steps = {}
                            for i in range(35):
                                if i ==0:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values
                                    start_+=8 #needed to begin the next index to keep up with proper dates
                                elif i<10:
                                    steps[i] = open_d10[f'{var_name}'][start_:start_+7,:,:].mean(dim=['step']).values #eight total possible values until last time step
                                    start_+= 8 #Need to add one because we don't want to re-index the same day
                                elif i == 10:
                                    #Need to take from the first file (time 00:00:00), and combine with d35 files
                                    s1 = (open_d35[f'{var_name}'][-1,:,:] + open_d35[f'{var_name}'][0,:,:] + \
                                        open_d35[f'{var_name}'][1,:,:] + open_d35[f'{var_name}'][2,:,:]) /4
                                    steps[i] = s1
                                    start_ = 3 #start count over, 4th file is the new date in d35 files
                                elif i <=34:
                                    steps[i] = open_d35[f'{var_name}'][start_:start_+4,:,:].mean(dim=['step']).values
                                    start_+=4
                            #Add to file
                            for step,lead_day in enumerate(steps.keys()):
                                template_GEFS_initial[:,idx,step,:,:] = steps[lead_day]           
                        #I can't save the S dimension, but it doesn't matter, we have the initialization date of the file already (it's the filename)
                        #Now all the data is into one file, so save as one file
           
#%%
            def julian_date(_date,template_GEFS_initial):
                #Return julian date for anomaly calculation
                a_date_in= template_GEFS_initial.shape[2]
                #get the start date
                a_start_date = pd.to_datetime(_date) 
    
                a_date_out=[]
                for a_i in range(a_date_in):
                    a_date_out.append((a_start_date + np.timedelta64(a_i,'D')).timetuple().tm_yday)
        
                return(a_date_out)

            #Add julian date
            julian_list = julian_date(_date,template_GEFS_initial)
            
            if 'apcp' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        pr = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
            
                    ),
                    attrs = dict(
                        Description = 'Total daily precipitation GEFSv12. Daily average already computed. All ensembles and Ls in one file'),
                )  
            elif 'dlwrf' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        dlwrf = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Downwelling longwave radiation GEFSv12. Daily average already computed. All ensembles and Ls in one file'),
                )  
                
           
            elif 'dswrf' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        dswrf = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Downwelling shortwave radiation GEFSv12. Daily average already computed. All ensembles and Ls in one file'),
                )  
            elif 'soil' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        RZSM = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
            
                    ),
                    attrs = dict(
                        Description = 'Volumetric soil moisture content at 4 levels: 0.0-0.1, 0.1-0.4, 0.4-1.0 and 1.-2. m depth \
    (fraction between wilting and saturation) GEFSv12. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'tmp' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        tmp = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Average temperature GEFSv12. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'ulwrf' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        ulwrf = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Longwave upwelling radiation. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'uswrf' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        uswrf = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Shortwave upwelling radiation. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'spfh' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        spfh = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Specific humidity. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'tmax' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        tasmax = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Maximum Temperature. Daily average already computed. All ensembles and Ls in one file')
                )
            elif 'tmin' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        tasmin = (['S','M','L','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_d10.X.values,
                        Y = open_d10.Y.values,
                        L = julian_list,
                        M = range(template_GEFS_initial.shape[1]),
                        S = np.atleast_1d(pd.to_datetime(_date)),
                    ),
                    attrs = dict(
                        Description = 'Minimum Temperature. Daily average already computed. All ensembles and Ls in one file')
                )
                
                
                
            # GEFS_out.assign_coords(S=out_date_create)
        
            out_name1 = 'out_1.nc4'
            GEFS_out.to_netcdf(path = f"{home_dir}/{var}/{out_name1}")
            #Now compress
            os.system(f'ncks -O -4 -L 1 {home_dir}/{var}/{out_name1} {save_out_dir_GEFS}/{final_out_name}4')
            #Remove the old files
            os.system(f'rm {home_dir}/{var}/{out_name1}')
            
#%%
if __name__ == '__main__':
    p=Pool(len(vars_to_process))
    p.map(merge_ensemble_members,vars_to_process)