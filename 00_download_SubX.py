#!/usr/bin/env python3

'''Create wget script for each SubX model to setup download of files in parallel 

Inputs: IRIDL NCAR username and password
Inputs -- you can change the model and variables that you want downloaded

Outputs: Shell script which contains download info. Files will download 20 at a time when run.

'''
import os
import datetime as dt
import numpy as np


username_IRIDL = ""
password_IRIDL = ""


home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX'
script_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts'

models = ['RSMAS']
sources = ['CCSM4']
vars = ['huss', 'dswrf','mrso','tas', 'uas', 'vas','tdps','pr','cape','tasmax','tasmin']

# get the dates for all SubX models (this is all the possible dates)
start_date = dt.date(1999, 1, 6)

#Testing end_date
#end_date = dt.date(1999, 1, 14)

#Actual end date
end_date = dt.date(2015, 12, 29)

#dates
dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

### GMAO GEOS_V2p1 model


#vars=['tas']

try:
    new_dir=(f'{home_dir}/{models[0]}')
    os.mkdir(new_dir)
except FileExistsError:
    print('Directory already created')

count=0    
output=[]
output.append('#!/bin/bash')
for m_i, model in enumerate(models):
    
    for d_i, date in enumerate(dates):
        
        for v_i, var in enumerate(vars):
        
            date_str = '{}-{}-{}'.format(str(date.year), str(date.month).rjust(2,'0'), str(date.day).rjust(2,'0'))

            command = f"wget -nc --user {username_IRIDL} --password {password_IRIDL} 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.{model}/.{sources[m_i]}/.hindcast/.{var}/S/(1200%20{str(date.day)}%20{date.strftime('%b')}%20{str(date.year)})/VALUES/data.nc' -O {home_dir}/{models[0]}/{var}_{model}_{date_str}.nc &"
            
            if count % 20 == 0:
                output.append('wait')
            
            count+=1
            output.append(command)
            
name = f'wget_{models[0]}'
np.savetxt(f'{script_dir}/{name}.txt',output, fmt="%s")



