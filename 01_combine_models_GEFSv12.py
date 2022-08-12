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

#for linux
home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX'
script_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts'
gefs_dir = f'{home_dir}/GEFSv12/raw_ensemble_files'

save_out_dir_GEFS = f'{home_dir}/GEFSv12/combined_ensembles_GEFS'
os.system(f'mkdir -p {save_out_dir_GEFS}')

mask_dir = f'{script_dir}/CONUS_mask'

vars_to_process= ['cape_sfc','dswrf_sfc', 'uflx_sfc','vflx_sfc', 'tmax_2m','tmin_2m', 'soilw_bgrnd', 'spfh_2m','pres_msl']

#Make all directories for save files (---neat list comprehension---)
[os.system(f'mkdir {save_out_dir_GEFS}/{i}') for i in vars_to_process]

#Test file

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
 
template_subx=xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/vas_GMAO_2015-02-09.nc4')
# #0-10 day
# template_gefs_10 = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GEFSv12/raw_ensemble_files/tmin_2m/d10_tmin_2m_2006091300_p02.nc4')
# #10-35 day
# template_gefs_35 = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GEFSv12/raw_ensemble_files/tmin_2m/d35_tmin_2m_2006091300_p02.nc4')

# desc='dswrf'

# var_t1 = xr.zeros_like(template_subx)

# #Get total length of lead days
# lead_all=len(template_gefs_10.step.values) + len(template_gefs_35.step.values)
# #now that we have an emtpy xarray, let's change the ensemble size and lead time
# model_num=11
# template_GEFS = np.empty(shape=(template_subx.S.shape[0],model_num,lead_all,template_gefs_10.Y.shape[0],template_gefs_10.X.shape[0]))
# template_GEFS.shape
#This template was created from the code above
template_GEFS_initial = np.empty(shape=(1,11,35,27,59))
lead_splices = ['d10','d35']

all_possible_ensemble_members = ['c00', 'p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09', 'p10']
#testing with one at a time
vars_to_process = ['pres_msl']
#Now begin adding data
for var in vars_to_process:
    # var=vars_to_process[0]
    print(f'Working on variable {var} to merge ensemble members.')
    os.chdir(f'{gefs_dir}/{var}')
    
    #Get the dates of the files
    for _date in dates:
        # _date=dates[0]

        path_out = f"{save_out_dir_GEFS}/{var}"
        out_date = f'{_date.year}-{_date.month:02}-{_date.day:02}'
        final_out_name= f"{var.split('_')[0]}_{out_date}.nc4"
        try:
            xr.open_dataset(f"{path_out}/{final_out_name}")
        except FileNotFoundError:
            template_GEFS_initial = np.empty(shape=(1,11,35,27,59))
            template_GEFS_initial[:,:,:,:,:] = np.nan
            lead_day=0 #keeps up with which index is correct in template_GEFS_initial (resets with each date)
            for splices in lead_splices:
                
                all_files_by_splice = sorted(glob(f'{splices}_{var}_{_date.year}{_date.month:02}{_date.day:02}00*'))
                #some files have doubles (rsync error or from HPC when converting)
                file_len = [len(i) for i in all_files_by_splice]
                mode = max(set(file_len), key=file_len.count)
                #Replace files
                all_files_by_splices = [i for i in all_files_by_splice if len(i) == mode]
                
                
                #If we have some realizations that are missing:
                if len(all_files_by_splice) < 11:
                    #Some ensembles are missing, split to get the name of ensemble members
                    avail_ensemble_members = [i.split('_')[-1].split('.')[0] for i in all_files_by_splice]
                    missing_members = list(set(all_possible_ensemble_members).difference(avail_ensemble_members))
                    
                    #Because we cannot go by strictly the number of the ensemble name (since p00 and c00 exist),
                    #we will just assign based on sorted order (which is all_possible_ensemble_members order)
                    for ensemble_number,ensem_zip in enumerate(zip(all_files_by_splice,all_possible_ensemble_members)):
                        #For missing ensemble members, we need to create an empty dataset to add to the file
                        # list(enumerate(zip(all_files_by_splice,all_possible_ensemble_members)))
                        if ensem_zip[1] not in missing_members:
                            #Now add data to file according to ensemble_number
                            open_f = xr.open_dataset(ensem_zip[0])
                            var_name = [i for i in list(open_f.keys()) if 'step' not in i][0] #get name of actual variable in file
                            #Now append to template_GEFS_initial
                            for step in range(open_f.step.shape[0]):
                                if var !='soilw_bgrnd':
                                    if splices == 'd10':
                                        template_GEFS_initial[:,ensemble_number,step,:,:] = open_f[f'{var_name}'][step,:,:]
                                    elif splices == 'd35':
                                        template_GEFS_initial[:,ensemble_number,11+step,:,:] = open_f[f'{var_name}'][step,:,:]
                                elif var == 'soilw_bgrnd':
                                    '''Important note, soilw_bgrnd has 4 levels. Need to take the sum of the first three levels through 0-1m)'''
                                    sum_soil_levels = ( open_f[f'{var_name}'][:,0,:,:] + \
                                        open_f[f'{var_name}'][:,1,:,:] + open_f[f'{var_name}'][:,2,:,:]).to_dataset()
                                    #Testing is fine
                                        # sum_soil_levels.soilw[0,0,0].values   
                                    # open_f[f'{var_name}'][0,0,0,0].values+\
                                    # open_f[f'{var_name}'][0,1,0,0].values+\
                                    # open_f[f'{var_name}'][0,2,0,0].values
                                    
                                    if splices == 'd10':
                                        template_GEFS_initial[:,ensemble_number,step,:,:] = sum_soil_levels[f'{var_name}'][step,:,:]
                                    elif splices == 'd35':
                                        template_GEFS_initial[:,ensemble_number,11+step,:,:] = sum_soil_levels[f'{var_name}'][step,:,:]
                                
                else:    
                    #If all possible files are there, then this is the easy code processing to add data to single file
                    for ensemble_number,file in enumerate(all_files_by_splice):
                        open_f = xr.open_dataset(file)
                        var_name = [i for i in list(open_f.keys()) if 'step' not in i][0] #get name of actual variable in file
                        #Now append to template_GEFS_initial
                        for step in range(open_f.step.shape[0]):
                            if var != 'soilw_bgrnd':
                                if splices == 'd10':
                                    #Files day 0-10 (11 days)
                                    template_GEFS_initial[:,ensemble_number,step,:,:] = open_f[f'{var_name}'][step,:,:]
                                elif splices == 'd35':
                                    #Files days 11-34 (or 35) (still not sure for all days)
                                    template_GEFS_initial[:,ensemble_number,11+step,:,:] = open_f[f'{var_name}'][step,:,:]
                            else:
                                '''Important note, soilw_bgrnd has 4 levels. Need to take the sum of the first three levels through 0-1m)'''
                                sum_soil_levels = ( open_f[f'{var_name}'][:,0,:,:] + \
                                    open_f[f'{var_name}'][:,1,:,:] + open_f[f'{var_name}'][:,2,:,:]).to_dataset()
                                if splices == 'd10':
                                    #Files day 0-10 (11 days)
                                    template_GEFS_initial[:,ensemble_number,step,:,:] = sum_soil_levels[f'{var_name}'][step,:,:]
                                elif splices == 'd35':
                                    #Files days 11-34 (or 35) (still not sure for all days)
                                    template_GEFS_initial[:,ensemble_number,11+step,:,:] = sum_soil_levels[f'{var_name}'][step,:,:]


            #I can't save the S dimension, but it doesn't matter, we have the initialization date of the file already (it's the filename)
            #Now all the data is into one file, so save as one file
            if 'tmin' in var:
                GEFS_out = xr.Dataset(
                    data_vars = dict(
                        tmin = (['S','model','lead','Y','X'], template_GEFS_initial[:,:,:,:,:]),
                    ),
                    coords = dict(
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
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
                      
                        X = open_f.X.values,
                        Y = open_f.Y.values,
                        lead = range(template_GEFS_initial.shape[2]),
                        model = range(template_GEFS_initial.shape[1]),
            
                    ),
                    attrs = dict(
                        Description = 'Meridional wind speed GEFSv12. Daily average already computed. All ensembles and leads in one file'),
                )   
    
        
        
            out_name1 = 'out_1.nc4'
            GEFS_out.to_netcdf(path = f"{gefs_dir}/{out_name1}")
            #Now compress
            os.system(f'ncks -O -4 -L 1 {gefs_dir}/{out_name1} {save_out_dir_GEFS}/{var}/{final_out_name}')
            #Remove the old files
            os.system(f'rm {gefs_dir}/{out_name1}')
#%%
