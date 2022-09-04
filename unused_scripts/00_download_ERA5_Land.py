#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ERA5-Land https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form

@author: kdl
"""



import cdsapi
import os
from multiprocessing import Pool
import numpy as np
import xarray as xr

n_processors = 31
os.chdir('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/ERA5_land/raw_data')

c = cdsapi.Client()



def multiprocess(day):
    month_range = []
    for i in range(1,13):
        month_range.append(f'{i:02}')
        
    var_list = ['2m_temperature', 'forecast_albedo', 'surface_net_solar_radiation',
    'total_precipitation', 'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2',
    'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4']

    for year in np.arange(2000,2023):
        for month in month_range:   
            for var in var_list:
                try:
                    xr.open_dataset(f'/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/ERA5_land/raw_data/{var}_{year}-{month}-{day}.nc')
                except FileNotFoundError:
                    c.retrieve(
                        'reanalysis-era5-land',
                        {
                            'variable': f'{var}',
                            'year': f'{year}',
                            'month': f'{month}',
                            'day': f'{day}',
                            'time': [
                                '00:00', '01:00', '02:00',
                                '03:00', '04:00', '05:00',
                                '06:00', '07:00', '08:00',
                                '09:00', '10:00', '11:00',
                                '12:00', '13:00', '14:00',
                                '15:00', '16:00', '17:00',
                                '18:00', '19:00', '20:00',
                                '21:00', '22:00', '23:00',
                            ],
                            'area': [
                                60, -150, 20,
                                -60,
                            ],
                            'format': 'netcdf',
                        },
                        f'/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/ERA5_land/raw_data/{var}_{year}-{month}-{day}.nc')
                    os.system('sleep 1')
                
    
day_range = []
for i in range(1,32):
    day_range.append(f'{i:02}')
    
if __name__ == "__main__":
    p=Pool(n_processors)
    p.map(multiprocess,day_range)