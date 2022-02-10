#!/ncar/usr/jupyterhub.hpc.ucar.edu/jupyterhub-20210616/bin/python

#Only get Soil Moisture from EMC GEFS hindcast

import os
import datetime as dt
import numpy as np


home_dir = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

# get the dates
start_date = dt.date(1999, 1, 6)
end_date = dt.date(2015, 12, 29)
models = ['EMC']
sources = ['GEFS']
vars = ['soilw1', 'soilw2']
dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

try:
    os.mkdir(f'{home_dir}/Data/SubX/')
    os.mkdir(f'{home_dir}/Data/SubX/{models[0]}')
except FileExistsError:
    print('Directory already created')
             

for m_i, model in enumerate(models):
    
    for d_i, date in enumerate(dates):
        
        for v_i, var in enumerate(vars):
        
            date_str = '{}-{}-{}'.format(str(date.year), str(date.month).rjust(2,'0'), str(date.day).rjust(2,'0'))

            command = "wget --user kdl0013@auburn.edu --password Mi$$ionary098 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.{}/.{}/.hindcast/.{}/S/({}%20{}%20{})/VALUE/data.nc' -O {}/Data/SubX/{}/{}_{}_{}.nc".format(sources[m_i], model, var, str(date.day), date.strftime('%b'), str(date.year), home_dir,model, var, model, date_str)
        
            os.system(command)
            
            
            
#%%
os.system(command)
