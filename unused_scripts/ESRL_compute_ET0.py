#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    EMC 
    

@author: kdl
"""

import xarray as xr
import numpy as np
from numba import njit,prange
import os
import datetime as dt
import pandas as pd



home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/EMC_test_data'
os.chdir(home_dir)

#Create date list for loop to compute ET0
start_date = dt.date(1999, 1, 6)
end_date = dt.date(2015, 12, 29)
dates = [str(start_date + dt.timedelta(days=d)) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

dates.index(DayofYear)


#To compute ET0, we need:
    #1.) Net Radiation (Rn) Mj m-2 day-1

def compute_ET0(DayofYear, date_list = dates):
    DayofYear = '1999-03-15'
    date_index = date_list.index(DayofYear)
    
    #Open 1 file to get proper dimensions of Model ----- Net Radiation
    Rn = xr.open_dataset(glob(f'*rad*{DayofYear}.nc')[0], engine='netcdf4', decode_cf=True)
    #Wind (u and v)
    windU = xr.open_dataset(glob(f'*uas*{DayofYear}.nc')[0], engine='netcdf4', decode_cf=True)
    windV = xr.open_dataset(glob(f'*vas*{DayofYear}.nc')[0], engine='netcdf4', decode_cf=True)
    #Temperature (units in Kelvin)
    Tmp = xr.open_dataset(glob(f'*tas*{DayofYear}.nc')[0], engine='netcdf4', decode_cf=True)
    #Specific Humidity (has 1 more dimensions P, S, M, L, Y, X) -- P is new dimension
    Sh = xr.open_dataset(glob(f'*huss*{DayofYear}.nc')[0], engine='netcdf4', decode_cf=True)
    
    Sh.huss[0,0,m,l,y,x].values
    
    


    #Square both u and v and then take the square root
    def compute_windSpeed(windU = windU, windV = windV):
        windSpeed = windU.copy()
        #Square u and v and take the square root to get windspeed
        windSpeed.uas[:,:,:,:,:] = xr.ufuncs.sqrt(windU.uas[:,:,:,:,:].values**2 + windV.vas[:,:,:,:,:].values**2)
        #Create a new name for file to not get confused
        windSpeed['windspeed'] = windSpeed.uas
        windSpeed = windSpeed.drop(['uas'])
        
        return (windSpeed)
    
    #Return Esat in Kpa 
    def compute_Esat(Tmp = Tmp):
        Esat = Tmp.copy()
        #Perform calculations
        Esat.tas[:,:,:,:,:] = (17.27 * (Esat.tas[:,:,:,:,:].values - 273.15)) / (237.3 + (Esat.tas[:,:,:,:,:].values  - 273.15))
        Esat.tas[:,:,:,:,:] = .611 * xr.ufuncs.exp(Esat.tas[:,:,:,:,:].values)
        Esat['esat'] = Esat.tas
        Esat = Esat.drop(['tas'])
        
        return(Esat)
        
        
        Esat.tas[:,:,:,:,:] =  xr.ufuncs.exp((17.27 * (Tmp.tas[:,:,:,:,:].values - 273.15) / 237.3 + (Tmp.tas[:,:,:,:,:].values  - 273.15)))
        
        Esat.tas[:,:,:,:,:] = Esat.tas[:,:,:,:,:] * .611
        
        exp(1) * 2
        
        Esat.tas[:,:,:,:,:].values 
        Esat2.tas[:,:,:,:,:].values      
        Tmp.tas[:,:,:,:,:].values
        
        xr.ufuncs.exp(Esat.tas[0,0,0,0,0].values)
        
        windSpeed.uas[:,:,:,:,:].values = sqrt(windSpeed.uas[:,:,:,:,:].values)
        
    #Number of realizations (M-dimension)
    Mshape = Rn.M.shape[0]
    #shape of forecast period (L-dimension) --- needed for dates
    Lshape = Rn.L.shape[0]
    #Y and X dimensions
    Yshape, Xshape = Rn.Y.shape[0], Rn.X.shape[0]
    
    #Compute ET0 for each M, L, Y, X (loop):
        #Save to a new netcdf file
        
    for m in range(Mshape):
        for l in range(Lshape):
            for y in range(Yshape):
                for x in range(Xshape):
                    
                    output = Rn.copy()
                    Rn.rad[0,m,l,y,x].values
    
    
 
    
    Rn[0,]
    #Constants
    
   
  
    
    Rn.info()
    Rn.L.info()
    Rn.S[:]
    pd.to_datetime(Rn.L[:])
    np.datetime64(Rn.L[:].values)
    Rn.variables
    
    
t1 = xr.open_dataset('tas_EMC_1999-01-06.nc')
t2 = xr.open_dataset('tas_EMC_1999-01-22.nc')

t1.tas[0,1,0,0,0].values
t2.tas[0,1,0,0,0].values
    
    
#%%Calculate ET0 

from glob import glob

for file in glob(f'{home_dir}/*.nc'):
    print(file)


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
