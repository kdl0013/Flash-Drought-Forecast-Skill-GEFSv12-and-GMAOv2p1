#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TODO: 
    Use the Evaporative Demand Drought Index definition and the Soil Moisture
    Percentile Drop definition for flash drought identification for all files:
        SubX - EDDI, RZSM
        EDDI - EDDI (EDDI)
        RZSM - RZSM (SMERGE)

@author: kdl
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
from multiprocessing import Pool


#open a single file for SubX reference ET and a single file from gridMET 
dir1 = 'main_dir'
num_processors = int('procs')

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

subX_dir = f'{dir1}/Data/SubX/GMAO'
subX_SM_out_dir = f'{dir1}/Data/SubX/GMAO/SM_converted_m3_m3' #conversion

output_dir = f'{gridMET_dir}/ETo_SubX_values' #Refernce ET output directory
smerge_in_dir = f'{dir1}/Data/SMERGE_SM/Raw_data' #raw files that have been boxed in CONUS
SM_SubX_out_dir = f'{dir1}/Data/SMERGE_SM/Raw_data/SM_SubX_values' #Smerge values overlayed on SubX grid
eddi_dir = f'{dir1}/Data/EDDI/convert_2_nc'
eddi_subX_dir = f'{dir1}/Data/EDDI/EDDI_SubX_values'

image_dir = f'{dir1}/Outputs/MSE_plots'

os.system(f'mkdir -p {image_dir}')
os.system(f'mkdir {SM_SubX_out_dir}')
os.system(f'mkdir {eddi_subX_dir}')
os.system(f'mkdir {output_dir}')

os.chdir(subX_dir) #Set directory for SubX
#Get date list for initialized files
date_list = sorted(glob('mrso*_*-*.nc'))
date_list = [i[-13:-3] for i in date_list]




