#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""



@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn
import netCDF4 as nc


dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
home_dir = f'{dir1}/Data/SubX/SubX_sm_grid_conversion'

mask_grd = f'{dir1}/Scripts/CONUS_mask/CONUS_mask.grd'

os.chdir(home_dir)
#%%
#read in file
gmao = np.loadtxt('GMAO_soil_moisture_meter_depth.txt')

gmao_df = pd.DataFrame(gmao, columns=['lon','lat','depth'])

gmao.reshape((360,180))

out_ = f'{home_dir}/gmao_depth.nc'
ds = nc.Dataset(out_,'w',format='NETCDF4')

lat = ds.createDimension('Y', len(np.unique(gmao[:,1])))
lon = ds.createDimension('X', len(np.unique(gmao[:,0])))

lats = ds.createVariable('Y', 'f4', ('Y',))
lons = ds.createVariable('X', 'f4', ('X',))
value = ds.createVariable('gmao_depth', 'f4', ('Y', 'X',))
value.units = 'meters (m)'

#Manually assign gmao, it didnt' work with the np.unique(values)
lats[:] = np.arange(-89.5,89.5+1,1)
lons[:] = np.arange(-179.5,179.5+1,1)

ds.close()

#now open with xarray to make it easier
new_open = xr.open_dataset(out_)
#%%
for lat in range(len(new_open.Y.values)):
    for lon in range(len(new_open.X.values)):
        
        out_Y, out_X = new_open.gmao_depth[lat,lon].coords.values()
        
        output = gmao_df[(gmao_df['lon'] == np.array(out_X)) & (gmao_df['lat'] == np.array(out_Y))]['depth']
        #replace values in file
        new_open.gmao_depth[lat,lon] = output.iloc[0]
#%%      
file_name = 'GMAO_depth_converted_from_dat.nc'
new_open.to_netcdf(f'{file_name}')

out_1 = 'GMAO_add_metdata.nc'
#Add metadata
os.system(f'ncatted -O -a units,X,c,c,"degrees_east" -a units,Y,c,c,"degrees_north" {file_name} {out_1}')

#Now convert to same as CONUS mask USDM grid
os.system(f'cdo -remapcon,{mask_grd} {out_1} GMAO_1x1_grid.nc ')

#complete
