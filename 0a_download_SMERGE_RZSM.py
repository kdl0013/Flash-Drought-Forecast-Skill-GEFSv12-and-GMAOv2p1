#!/usr/bin/env python3

'''Downloading data for SMERGE soil moisture
Inputs will come from shell script 

Step 1, download links list from website 
https://disc.gsfc.nasa.gov/datasets/SMERGE_RZSM0_40CM_2.0/summary


Step 2, run script to download all files in parallel. Start year for download 
is 1999.

Notes: 
    This script will know if the files are already downloaded so it will not
    download these files again
'''

import os
import numpy as np
import requests
from multiprocessing import Pool

home_dir = 'main_dir'
num_processors = int('procs')
#%%
# #for testing (don't run)
home_dir='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
num_processors=12

#GES DISC link list data 
link_list_name='subset_SMERGE_RZSM0_40CM_2.0_20220120_141903.txt'

#Make directories
#For all data
try:
    os.mkdir(os.path.join(home_dir,'Data'))
    os.mkdir(os.path.join(home_dir,'Data','SMERGE_SM'))
    os.mkdir(os.path.join(home_dir,'Data','SMERGE_SM','Raw_data'))
except FileExistsError:
    print('Created directory folders')


#Read link list file
smerge_data = np.loadtxt(os.path.join(home_dir,'Data','SMERGE_SM',link_list_name),skiprows=1,dtype='str')

'''Instructions from GES DISC website https://disc.gsfc.nasa.gov/data-access#python-requests

  # Set the URL string to point to a specific data URL. Some generic examples are:
   #   https://servername/data/path/file
   #   https://servername/opendap/path/file[.format[?subset]]
   #   https://servername/daac-bin/OTF/HTTP_services.cgi?KEYWORD=value[&KEYWORD=value]
   URL = 'your_URL_string_goes_here'
   
   # Set the FILENAME string to the data file name, the LABEL keyword value, or any customized name. 
   FILENAME = 'your_filename_string_goes_here'
   
   import requests
   result = requests.get(URL)
   try:
      result.raise_for_status()
      f = open(FILENAME,'wb')
      f.write(result.content)
      f.close()
      print('contents of URL written to '+FILENAME)
   except:
      print('requests.get() returned an error code '+str(result.status_code))
 '''


#save data to correct path, function will check if file is already downloaded
os.chdir(os.path.join(home_dir,'Data','SMERGE_SM','Raw_data'))

already_downloaded_files = os.listdir()

def download_SMERGE(link_list_URL):
    if (link_list_URL[-12:] in already_downloaded_files):
        print(f'{link_list_URL[-12:]} already downloaded.')
    else:
        if (int(link_list_URL[-12:-8]) < 1999) or (int(link_list_URL[-12:-8]) > 2016):
            print(f'{link_list_URL[-12:]} year not needed.')
        else:
            result = requests.get(link_list_URL)
            try:
               result.raise_for_status()
               f = open(link_list_URL[-12:],'wb')
               f.write(result.content)
               f.close()
               print('contents of URL written to '+link_list_URL[-12:])
            except:
               print('requests.get() returned an error code '+str(result.status_code))
    


print('Downloading SMERGEv2.0 soil moisture data from GEOS DISC')

if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(download_SMERGE,smerge_data)


# #No multiprocessing (much slower)
# for URL in smerge_data:
#     result = requests.get(URL)
#     try:
#         result.raise_for_status()
#         f = open(URL[-12:],'wb')
#         f.write(result.content)
#         f.close()
#         print('contents of URL written to '+URL[-12:])
#     except:
#         print('requests.get() returned an error code '+str(result.status_code))

