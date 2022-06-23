#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
I noticed that some files were blank after ETo and RZSM anomaly processing e.g., 2004-06-24.

But ETo and RZSM anomaly for 2001-06-24 is fine.
@author: kdl
"""

import xarray as xr
import numpy as np
from glob import glob
import os

dir1 = 'main_dir'
model_NAM = 'model_name'

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM = 'GMAO'

home_dir = f'{dir1}/Data/SubX/{model_NAM}'
os.chdir(home_dir)

subX_dir = f'{dir1}/Data/SubX/{model_NAM}/' #where subX model data lies 
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask file
    

'''

USDM-West_mask:
---------------

The U.S. Drought Monitor (USDM) Regions
 -- For the NCA-LDAS domain

Data description:

 -- This dataset includes CONUS-wide (minus Alaska and Hawaii for now)
    U.S. Drought Monitor regions, delimited by state-political boundaries.
    These regions were defined by the U.S. Drought Monitor.

Legend Information:

 -- Each region has been assigned an integer-based index value here.
    The corresponding region and integer values include:
    1:  West region
    2:  Midwest region
    3:  HighPlains region
    4:  South region
    5:  Southeast region
    6:  Northeast region

 -- Note: There are two separate USDM masks - one which should be used
     for the HighPlains region and one that should be used for the West
     region.  The reason for this is that the states of Colorado and
     Wyoming are in BOTH of these regions, as defined by the USDM.
     Thus, any analysis of the "HighPlains" region while using the West
     mask will be INCORRECT.  Same for an analysis of the "West" region
     while using the HighPlains mask.  The spatial extents for all other
     regions are identical between the two different masks.


For more references and information, please visit:
 -- https://www.drought.gov/
 -- https://droughtmonitor.unl.edu/


'''




#Get date list
test_var='ETo'

def return_date_list():
    date_list = []
    for file in sorted(glob(f'{home_dir}/{test_var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list() 

#%%
'''Step 1
re-open anomaly files (ETo and RZSM) and then add the mask for np.nan based on HP_conus_mask
'''

var='RZSM'

for date in init_date_list:
    file_list = f'{var}_anomaly_{date}.nc4'
    a = xr.open_dataset(file_list)
    a.close()
    
    for Y in range(a.Y.shape[0]):
        for X in range(a.X.shape[0]):
            if HP_conus_mask.High_Plains[0,Y,X].values in range(1,7):
                a[f'{var}_anom'][:,:,:,Y,X] = np.nan
                
    a.to_netcdf('test.nc4')



#%%
'''Now we want to open up each ETo- and RZSM anomaly file (also EDDI), to see
which ones have missing data'''
var = 'RZSM'
file_list = f'{var}_anomaly_*.nc4'

sorted(glob(file_list))
#test with single file to see 
#RZSM_anomaly
missing_anomaly_RZSM= {} #subx missing anomaly data files
missing_original_RZSM = [] #mrso missing data files
missing_originial_RZSM_m3_m3 = [] #missing converted mrso m3/m3 files
missing_SMERGE_to_SubX_files = [] #when we convert smerge format to subx format for ACC comparison


for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    open_f.close()
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    #anomaly. I inspected a file that was filled correctly and it had 528 missing value
    if np.count_nonzero(open_f.RZSM_anom[0,:,:,:,:].values < 0) < 528:
        missing_anomaly_RZSM[file] = np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values)) 
        
    open_mrso = xr.open_dataset(f'mrso_{model_NAM}_{file[-14:-4]}.nc') 
    if np.count_nonzero(np.isnan(open_mrso.mrso[0,:,:,:,:].values)) == 286740:
        missing_original_RZSM.append(file)

    open_m3_m3 = xr.open_dataset(f'SM_converted_m3_m3/SM_SubX_m3_m3_{file[-14:-4]}.nc4') 
    if np.count_nonzero(np.isnan(open_m3_m3.SM_SubX_m3_m3_value[0,:,:,:,:].values)) == 286740:
        missing_originial_RZSM_m3_m3.append(file)
    
    #SMERGE only
    open_SMERGE_subx = xr.open_dataset(f'{dir1}/Data/SMERGE_SM/SM_SubX_values/SM_SubX_{file[-14:-4]}.nc') 
    if np.count_nonzero(np.isnan(open_m3_m3.SM_SubX_m3_m3_value[:,:,:,:,:].values)) < 80640:
        missing_SMERGE_to_SubX_files.append(file)        


print(f'Missing data from {len(missing_anomaly_RZSM)} RZSM anomaly files')
print(f'Missing data from {len(missing_original_RZSM)} RZSM original (subx mrso file) files')
print(f'Missing data from {len(missing_originial_RZSM_m3_m3)} RZSM m3/m3 subx files')
print(f'Missing data from {len(missing_SMERGE_to_SubX_files)} RZSM m3/m3 subx files')

#%%
var = 'ETo'

missing_anomaly_ETo = []
missing_original_ETo = []

if var == 'EDDI':
    file_list = f'{var}_*.nc4'
else:
    file_list = f'{var}_anomaly_*.nc4'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    if np.count_nonzero(np.isnan(open_f.ETo_anom[0,:,:,:,:].values)) == 138240:
        missing_anomaly_ETo.append(file)
        
    open_ETo = xr.open_dataset(f'{var}_{file[-14:-4]}.nc4') 
    if np.count_nonzero(np.isnan(open_ETo.ETo[0,:,:,:,:].values)) == 286740:
        missing_original_ETo.append(file)

print(f'Missing data from {len(missing_anomaly_ETo)} ETo anomaly files')
print(f'Missing data from {len(missing_original_ETo)} ETo original (subx file)')

# missing_anomaly_ETo
# missing_original_ETo
#%%           
#No anomalies for EDDI
var = 'EDDI'

missing_original_EDDI = []

if var == 'EDDI':
    file_list = f'{var}_*.nc4'
else:
    file_list = f'{var}_anomaly_*.nc4'

for file in sorted(glob(file_list)):
    open_f = xr.open_dataset(file)
    '''after further inspection, if a file has all nan values then it has a total
    of 286740 after function np.count_nonzero(np.isnan(open_f.RZSM_anom[0,:,:,:,:].values));
    therefore, this is a missing data file.
    '''
    if np.count_nonzero(np.isnan(open_f.EDDI[0,:,:,:,:].values)) == 286740:
        missing_original_EDDI.append({f'{file}':np.count_nonzero(np.isnan(open_f.EDDI[0,:,:,:,:].values))})
        
print(f'Missing data from {len(missing_original_EDDI)} EDDI subx files')

# missing_original_EDDI


