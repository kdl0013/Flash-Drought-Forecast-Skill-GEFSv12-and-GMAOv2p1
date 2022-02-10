#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing my capabilities with xarray and playing with more dimensions

@author: kdl
"""

import xarray as xr
import numpy as np
from numba import njit,prange
import os
import datetime

file = xr.open_dataset('/home/kdl/Downloads/wtas_USA2_1999-01-06.nc')
file.info()


os.chdir('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/EMC_test_data/append_test')
print(os.listdir())
#%%Open up single file, soilw1


    

#Create a soil class that works with only soil files with SubX data
class soil_EMC():
    def __init__(self):
        self = self
    
    #Read file, create a new variable (soil_total) which adds together soilw1 & 2
    def read_file(self):
        file = xr.open_dataset(self, decode_cf = True)
        file = file.assign(soil_total = file.soilw1 + file.soilw2)
        file_start_date = pd.to_datetime(file.S.values).date
        
        return(file, file_start_date)

    def RZSM_percentile(glob_directory):    

        


soil_nc, file_start = soil_EMC.read_file('soilw2_EMC_2004-11-13.nc_compressed.nc')
soil_nc
#%%


#To construct a climatology, I need to make a percentile for each:
    #Grid point:
        #Lead time:
            #Model


soil_nc.soil_total[0,1,1,50,90].values
soil_nc.soil_total[0,1,1,5:,:].shape
pd.to_datetime(soil_nc.L[:])



a = pd.to_datetime(soil_nc.S.values).date
a


soil_nc.L.day
datetime.


a = xr.decode_cf(soil_nc)
a.L[:]

np.datetime64(soil_nc.L[:].values)

soil1 = xr. open_dataset('soilw1_EMC_2001-08-22.nc_compressed.nc')

soil1.info()
soil1.dimensions()

test
