#!/usr/bin/env python3

'''Create wget script for each SubX model to setup download of files in parallel.

This can only be done with GMAO GEOS. ESRL FIM1r1p1 & RSMAS CCSM4 do not have the 
correct near surface variables.  EMC GEFSv12 requires downloading with a different script. 

Inputs: IRIDL NCAR username and password
Inputs -- you can change the model and variables that you want downloaded

Outputs: Shell script which contains download info. Files will download 20 at a time when run.



!!!FOR GMAO, only every 5 days are initiated beginnning January 10, 1999
'''
import os
import datetime as dt
import numpy as np
from glob import glob


username_IRIDL = ""
password_IRIDL = ""


home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX'
script_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Scripts'
model = 'GMAO'

#Works for GMAO GEOSo
models = ['GMAO']
sources = ['GEOS_V2p1']

model_dir = f'{home_dir}/{model}'

vars = ['huss', 'dswrf','mrso','tas', 'uas', 'vas','tdps','pr','cape','tasmax','tasmin']

# get the dates for all SubX models (this is all the possible dates)
#GMAO starts at 1999, 1, 10

#RSMAS is only initialized on Sundays
start_date = dt.date(1999, 1, 6)

#Testing end_date
#end_date = dt.date(1999, 1, 14)

#Actual end date
end_date = dt.date(2015, 12, 29)

#dates
dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

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
    
    if model == 'GMAO':
        #dates
        dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]
        variables = ['huss', 'dswrf','mrso','tas', 'uas', 'vas','tdps','pr','cape','tasmax','tasmin']

    # elif model == 'RSMAS':
    #     dates = [i for i in dates if i.weekday()==6] #RSMAS only initialized on Sundays
    #     variables = ['huss', 'dswrf','mrso','tas', 'uas', 'vas','tdps','pr','cape','tasmax','tasmin']

    for d_i, date in enumerate(dates):
        
        for v_i, var in enumerate(variables):
            # if count == 60:
            #     break
            date_str = '{}-{}-{}'.format(str(date.year), str(date.month).rjust(2,'0'), str(date.day).rjust(2,'0'))
            
            #GMAO
            # command = f"wget -nc --user {username_IRIDL} --password {password_IRIDL} 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.{model}/.{sources[m_i]}/.hindcast/.{var}/S/(1200%20{str(date.day)}%20{date.strftime('%b')}%20{str(date.year)})/VALUES/data.nc' -O {home_dir}/{models[0]}/{var}_{model}_{date_str}.nc &"
            
            #RSMAS
            command = f"wget -nc --user {username_IRIDL} --password {password_IRIDL} 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.{model}/.{sources[m_i]}/.hindcast/.{var}/S/(1200%20{str(date.day)}%20{date.strftime('%b')}%20{str(date.year)})/VALUES/data.nc' -O {home_dir}/{models[0]}/{var}_{model}_{date_str}.nc4 &"

            if count % 20 == 0:
                output.append('wait')
            
            count+=1
            output.append(command)
            
   
name = f'wget_{models[0]}'
np.savetxt(f'{script_dir}/{name}.txt',output, fmt="%s")

os.system(f'mv {script_dir}/{name}.txt {script_dir}/{name}.sh')

#%%
#make individual dates 
#first find missing dates for dswrf (five were missing some data)
model_dir = f'{home_dir}/{model}'
os.chdir(model_dir)
all_dates = [date_[-14:-4] for date_ in sorted(glob('mrso*.nc4'))]
dswrf_dates = [date_[-14:-4] for date_ in sorted(glob('dsw*.nc4'))]

missing_dates=set(all_dates).difference(dswrf_dates)
    
def file_list():
    return([i for i in glob('*.nc4')])
file_list()
