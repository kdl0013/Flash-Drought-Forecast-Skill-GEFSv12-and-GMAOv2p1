#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix NRL issue with dates that was not fixed with the 01e_remove_S_dimension script


@author: kdl
"""

import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob
import bottleneck as bn
from multiprocessing import Pool
from numpy import inf
import sys
import datetime

dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
os.chdir(f'{dir1}/Data/SubX/NRL')



var_list=["actual_vapor_pressure" ,"srad",  "windspeed", "tas"]
def change_S_dim_again(var):
# for var in var_list:
    print(f'Working on var {var}')
    for f in sorted(glob(f'{var}*.nc4')):
        try:
            a=xr.open_dataset(f)
        except OSError:
            os.system(f'rm {f}')
        doy = []
        for d in range(len(a.L.values)):
            d1=pd.to_datetime(a.S.values[0]) + np.timedelta64(a.L.values[d],'D')
            doy.append(pd.to_datetime(d1).date().timetuple().tm_yday)
            
        a['L'] = doy
        bb=pd.to_datetime(a.S.values[0]).date()
        out_test = f'{bb.year}-{bb.month:02}-{bb.day:02}'
        
        a['S'] = np.atleast_1d(out_test)
        a.to_netcdf(path='a_test.nc4')
        a.close()
        os.system(f'mv a_test.nc4 {f}')
#%%
if __name__ == '__main__':
    p=Pool(len(var_list))
    p.map(change_S_dim_again,np.array(var_list))