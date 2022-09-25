#!/usr/bin/env python3

'''Because all GEFSv12 models are in seperate files for each model, I need to save
on space. So save grib2 files as a netcdf, convert data to CONUS grid cell,
compress, and re-send back to home computer.

Run this on HPC since that is where data is housed

https://noaa-gefs-retrospective.s3.amazonaws.com/Description_of_reforecast_data.pdf

'''
import os
import datetime as dt
import numpy as np
import xarray as xr
from glob import glob
from multiprocessing import Pool

# #for Easley cluster
scratch_dir = '/scratch/kdl0013' #where all directories are stored
home_dir = '/home/kdl0013/GEFSv12'
mask_dir = '/home/kdl0013'

#Test file
# file_o = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GEFSv12/raw_ensemble_files/pres$

'''Compress file size to send back to home computer'''
def cluster_compression_for_GEFS_files(year):
    var = 'ulwrf_sfc'
    os.chdir(f'{scratch_dir}/{var}')
    for file in sorted(glob('*{year}*.grib2'),reverse=True):
        file_o = file

        try:
            save_dir = f'{home_dir}/{var}'
            xr.open_dataset(f'{save_dir}/{file_o[:-6]}.nc4')
        except FileNotFoundError :
            #Some ensemble member files are missing (I can phsyically see this in the directory)
                #Still need to figure out what to do with these files once they are read into this$
#TODO:

            #Open file
            # var='cape_sfc'
            # file_o = os.listdir()[20]
            #Because there are different files, we need to use the proper opening format for each $
            '''Either cf or pf, depending on the file name'''

            var = file_o.split('_')[1] + '_' + file_o.split('_')[2]

            if file_o.split('_')[-1].split('.')[0] == 'c00':
                try:
                    grib_o = xr.open_dataset(file_o,filter_by_keys={'dataType': 'cf'}) #This works
                    #save to netcdf to compute daily mean for each variable
                    file_nc_name = f"{file_o[:-6]}.nc4"
                    var_location = f'{scratch_dir}/{var}'
                    grib_o.to_netcdf(f"{var_location}/{file_nc_name}")
                    #Now use cdo operators on the file. Daily mean, remap to CONUS grid, select only CONUS$
                    #remove old nc file (to save on storage space)

                    os.system(f'cdo -remapcon,{mask_dir}/CONUS_mask.grd  {var_location}/{file_nc_name} {save_dir}/{file_nc_name}')
                    os.system(f"rm {var_location}/{file_nc_name}")
                    
                except EOFError:
                    pass #error caused by realization having no data
                except ValueError:
                    pass #error caused by realization having no data
            else:
                try:
                    grib_o = xr.open_dataset(file_o,filter_by_keys={'dataType': 'pf'}) #This works$
                    file_nc_name = f"{file_o[:-6]}.nc4"
                    var_location = f'{scratch_dir}/{var}'
                    grib_o.to_netcdf(f"{var_location}/{file_nc_name}")
                    #Now use cdo operators on the file. Daily mean, remap to CONUS grid, select only CONUS$
                    #remove old nc file (to save on storage space)

                    os.system(f'cdo -remapcon,{mask_dir}/CONUS_mask.grd {var_location}/{file_nc_name} {save_dir}/{file_nc_name}')
                    os.system(f"rm {var_location}/{file_nc_name}")
                except EOFError:
                    pass
                except ValueError:
                    pass #error caused by realization having no data

# #Compress once more
# final_outname = f"{file_o.split('/')[-1][:-6]}.nc4"
# os.system(f'ncks -4 -L 1 {var_location}/{remap_nc_name} {var_location}/{final_outname}')
# #Remove the old files
# os.system(f'rm {var_location}/{file_nc_name} && rm {var_location}/{remap_nc_name}')
# #Move final output into final save location
# os.system(f'cp {var_location}/{final_outname} {save_dir}/{final_outname}')


#vars_to_process= ['uswrf_sfc']


if __name__ == '__main__':
    p = Pool(5)
    p.map(cluster_compression_for_GEFS_files,np.arange(2000,2023))


