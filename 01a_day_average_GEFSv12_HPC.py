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

dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'

# #for linux
home_dir = f'{dir1}/Data/SubX/GEFSv12/raw_ensemble_files/'
out_dir = f'{dir1}/Data/SubX/GEFSv12/combined_ensembles_GEFS/'

vars_to_process= ['cape_sfc','dswrf_sfc', 'uflx_sfc','vflx_sfc','soilw_bgrnd', 'tmax_2m','tmin_2m', 'spfh_2m','pres_msl']

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

for var in vars_to_process:
    # var=vars_to_process[0]
    print(f'Working on variable {var} to merge ensemble members.')
    os.chdir(f'{home_dir}/{var}')
    save_out_dir_GEFS = f'{out_dir}/{var}'
    
    #Get the dates of the files
    for _date in dates:
        # _date=dates[0]

        path_out = f"{save_out_dir_GEFS}"
        #The date of the file is when it was intialized
        out_date_create = pd.to_datetime(_date) 
        out_date = f'{out_date_create.year}-{out_date_create.month:02}-{out_date_create.day:02}'
                
        
        final_out_name= f"{var.split('_')[0]}_{out_date}.nc4"
        try:
            xr.open_dataset(f"{path_out}/{final_out_name}")
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
               
                
            #If we have some realizations that are missing:
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

            
            elif (len(all_files_d10) == 0 )and (len(all_files_d35) == 0):
                pass
            
            elif (len(all_files_d10) < 11) and (len(all_files_d35) < 11):
 
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
#%%             
            #I can't save the S dimension, but it doesn't matter, we have the initialization date of the file already (it's the filename)
            #Now all the data is into one file, so save as one file
            
            #Only look at complete files first,
            if True:
                if 'tmin' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            tmin = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Tmin GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                elif 'tmax' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            tmax = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Tmax GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                    
                elif 'cape' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            cape = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Convective active potential energey GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                elif 'dswrf' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            dswrf = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Downwelling shortwave radiation GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                elif 'pres' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            SLP = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Mean sea level pressure GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                elif 'soil' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            RZSM = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Volumetric soil moisture content at 4 levels: 0.0-0.1, 0.1-0.4, 0.4-1.0 and 1.-2. m depth \
        (fraction between wilting and saturation) GEFSv12. Daily average already computed. All ensembles and leads in one file')
                    )
                elif 'spfh' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            humidity = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Specific humidity GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )  
                elif 'uflx' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            u_wind = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Zonal wind speed GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )           
                elif 'vflx' in var:
                    GEFS_out = xr.Dataset(
                        data_vars = dict(
                            v_wind = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                        ),
                        coords = dict(
                          
                            X = open_d10.X.values,
                            Y = open_d10.Y.values,
                            lead = range(template_GEFS_initial.shape[2]),
                            model = range(template_GEFS_initial.shape[1]),
                
                        ),
                        attrs = dict(
                            Description = 'Meridional wind speed GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                    )   
        
                # GEFS_out.assign_coords(S=out_date_create)
            
                out_name1 = 'out_1.nc4'
                GEFS_out.to_netcdf(path = f"{home_dir}/{var}/{out_name1}")
                #Now compress
                os.system(f'ncks -O -4 -L 1 {home_dir}/{var}/{out_name1} {save_out_dir_GEFS}/{final_out_name}')
                #Remove the old files
                os.system(f'rm {home_dir}/{var}/{out_name1}')
#%%
