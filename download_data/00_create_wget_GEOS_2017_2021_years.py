#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Location of GMAO GEOS-5 hindcasts https://gmao.gsfc.nasa.gov/gmaoftp/gmaofcst/subx/GEOS_S2S_V2.1_fcst/IRI/

File name order: cape_GMAOGEOS_01apr2018_00z_d01_d45_m01.nc 
--- Only get every 5 days. Starts on 25jul2017.
---Need to get all 4 models

@author: kdl
"""

import xarray as xr
import numpy as np
import pandas as pd
import datetime
import os

dir1='/home/kdl/Insync/OneDrive/NRT_CPC_Internship/'
# dir1 = 'main_dir'
out_dir = f'{dir1}/Data/SubX/GMAO/preprocess_2017_2021'
vars_to_download = ['dswrf_sfc', 'uas_10m','vas_10m', 'tasmax_2m','tasmin_2m', 'mrso', 'spfh_2m','cape_sfc','tdps_2m','huss_850']

#
start_date = pd.to_datetime('07-25-2017')
end_date = pd.to_datetime('12-31-2021')

url0 = 'https://gmao.gsfc.nasa.gov/gmaoftp/gmaofcst/subx/GEOS_S2S_V2.1_fcst/IRI/'
output=[]
count=0
for i in range(400): #random rannge to get all dates between 2017-2021
    # start_date = start_date
    
    if start_date < end_date:

    
        def convert_date(_date):
            day = f'{start_date.day:02}'
            year = start_date.year
            month = _date.strftime('%B').lower()[0:3]
            out_date = f'{day}{month}{year}'
            return(out_date)
        
        out_date = convert_date(start_date)
        
        def output_file_name(var,out_date,model):
            varname = var.split('_')[0] + '_GMAO_'
            date_name = f'{pd.to_datetime(out_date).year}-{pd.to_datetime(out_date).month:02}-{pd.to_datetime(out_date).day:02}_m{model}.nc4'
            out_name = f'{varname}{date_name}'
            return(out_name)
        
        #add to url with variable
        for var in vars_to_download:
            for model in [1,2,3,4]:
                if count%30 == 0:
                    output.append('wait')
                    
                url1 = f'{url0}{var}_GMAOGEOS_{out_date}_00z_d01_d45_m{model:02}.nc'
                out_name = output_file_name(var,out_date,model)
                # output.append(f'wget -nc -O {url1} {out_dir}/{out_name} &') #doesn't work
                output.append(f'wget -nc {url1} &')
                count+=1
    #Now only get data for every 5 days
    start_date = start_date + np.timedelta64(5,'D')
        
        
                    
np.savetxt(f'{dir1}/Scripts/wget_GMAO_2021.txt',output, fmt="%s")

os.system(f"cat {dir1}/Scripts/wget_GMAO_2021.txt > {dir1}/Scripts/wget_GMAO_2021.sh")
os.system(f"rm {dir1}/Scripts/wget_GMAO_2021.txt")

        
        
        
