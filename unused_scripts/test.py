# ! /ncar/usr/jupyterhub.hpc.ucar.edu/jupyterhub-20210616/bin/python

import os
import datetime as dt


main_dir='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'

while True:
    answer = input(f'Is this the correct main directory as previous scripts  {main_dir}  ? YES/NO ').upper()
    if answer=='YES':
        break
    if answer=='NO':
        main_dir=input('What is the correct directory?')
    if answer != 'YES' and answer !='NO':
        continue

        
# get the dates
start_date = dt.date(1999, 1, 6)
end_date = dt.date(1999, 1, 8)
#end_date = dt.date(2015, 12, 29)
models = ['EMC']
sources = ['GEFS']
vars = ['soilw1', 'soilw2']
dates = [start_date + dt.timedelta(days=d) for d in range(0, end_date.toordinal() - start_date.toordinal() + 1)]

try:
    os.mkdir(os.path.join(main_dir,f'Data/SubX/{models[0]}_{sources[0]}'))
except FileExistsError:
    print('Directory already created')

for m_i, model in enumerate(models):
    
    for d_i, date in enumerate(dates):
        
        for v_i, var in enumerate(vars):
        
            date_str = '{}-{}-{}'.format(str(date.year), str(date.month).rjust(2,'0'), str(date.day).rjust(2,'0'))

            command = f"wget --user kdl0013@auburn.edu --password Mi$$ionary098 'http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/.{model}/.{sources[m_i]}/.hindcast/.{var}/S/(1200%20{str(date.day)}%20{date.strftime('%b')}%20{str(date.year)})/VALUE/data.nc' -O /home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/{model}_{sources[0]}/{var}_{model}_{date_str}.nc"
        
            os.system(command)
        
print(command)