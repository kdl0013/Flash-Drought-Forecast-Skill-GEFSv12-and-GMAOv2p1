#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Download Weekly EDDI Data from NOAA.gov EDDI website
Function will download all weekly data and save to the directory that is 
requested when running the script.

Inputs: Directory name
Outputs: Weekly EDDI anomaly files for North America

Best if this script is run through in Linux distribution

@author: Kyle Lesinger
"""

from calendar import monthrange
import os
import numpy as np


home_dir = 'main_dir'
script_dir = f'{home_dir}/Scripts'
#Create a file for multiproc

os.chdir(script_dir)

# home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

date_list = []
date_list.append('#!/bin/bash')
#Date ranges                    
years = range(1999, 2017)
months = range(1, 13) 
#Loop for file download where EDDI data is saved
url0 = 'https://downloads.psl.noaa.gov/Projects/EDDI/CONUS_archive/data/'
   
for y in years:
    count = 0
    url1 = url0 + str(y) + "/"
    for m in months:
        ym = str(y)+ str("{:02d}".format(m))
        for d in range(1,monthrange(y,m)[1]+1): #with monthrange you have the days in a month
            url2 = url1 + "EDDI_ETrs_01wk_" + ym + str("{:02}".format(d)) +".asc"
            #Make wget list
            wget = f'wget -nc -c -nd {url2} &'
            
            date_list.append(wget)
            
            if count % 30 ==0:
                date_list.append('wait')
            
            count +=1

np.savetxt(os.path.join(script_dir,'0d_download_EDDI_wget.sh'), date_list, fmt = '%s')                                       
