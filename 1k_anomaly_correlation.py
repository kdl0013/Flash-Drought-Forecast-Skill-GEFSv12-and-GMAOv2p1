#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TODO: Find the pearson_r correlation initially for:
    1.) The entire CONUS
        a.) for each model seperately
        b.) for each weekly lead time
        c.) plot

    2.) For each region individually
        a.) model, weekly lead time, season
@author: kdl
"""

import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.style.use('seaborn-talk')

import xarray as xr
import numpy as np
from numba import njit, prange
import numpy.ma as ma
import os
import datetime as dt
import pandas as pd
from glob import glob
from scipy.stats import rankdata
import xskillscore as xs
import sys
import dask
from climpred import HindcastEnsemble
import climpred
import cartopy.crs as ccrs
import cartopy as cart
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib as mpl

import warnings
warnings.filterwarnings("ignore")

print(f"PYTHON: {sys.version}")  # PYTHON: 3.8.1 | packaged by conda-forge | (default, Jan 29 2020, 15:06:10) [Clang 9.0.1 ]
print(f" xarray {xr.__version__}")  # xarray 0.14.1
print(f" numpy {np.__version__}")  # numpy 1.17.3
print(f" matplotlib {mpl.__version__}")  # matplotlib 3.1.2
print(f" cartopy {cartopy.__version__}")  # cartopy 0.17.0

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# model_NAM1 = 'GMAO'

os.chdir(f'{dir1}/Data/SubX/{model_NAM1}')


#TODO change later
# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
# mod = 'GMAO'
dir1 = 'main_dir'
model_NAM1 = 'model_name'
subX_dir = f'{dir1}/Data/SubX/{model_NAM1}/' #where subX model data lies 
HP_conus_path = f'{dir1}/Data/CONUS_mask/High_Plains_mask.nc'
West_conus_path = f'{dir1}/Data/CONUS_mask/West_mask.nc'

#observations
obs_eto_path = f'{dir1}/Data/gridMET/ETo_SubX_values/'
obs_rzsm_path = f'{dir1}/Data/SMERGE_SM/SM_SubX_values/'

output_nc_dir = f'{subX_dir}/skill_assessments/' #for skill assessment .nc files
output_image_dir = f'{dir1}/Outputs/pearson_correlation/'
output_season_dir = f'{dir1}/Outputs/pearson_correlation/seasonal_skill/'

'''Weekly pearson R skill score
Read all files in at once with xr.open_mfdataset, then calculate the skill
based only on weekly lead time. '''

var = 'RZSM'
# sub_name='vas'
def all_skill_by_lead (var, obs_eto_path, obs_rzsm_path):
    
    try:
        var_yearly = xr.open_dataset(f'{output_nc_dir}/interannual_pearson_skill.nc')
        return(var_yearly)

    except FileNotFoundError:
            
        if var == 'ETo':
            obs_name = 'ETo_SubX_*.nc'
            sub_name = 'ETo_anom'
            obs_path = obs_eto_path
        elif var == 'RZSM':
            obs_name = 'SM_SubX_*.nc'
            sub_name = 'RZSM_anom'
            obs_path= obs_rzsm_path
        
        subx_files = xr.open_mfdataset(f'{var}_anomaly*.nc', concat_dim=['S'], combine='nested')
        obs_files = xr.open_mfdataset(f'{obs_path}/{obs_name}', concat_dim = ['S'], combine = 'nested')
        
        #Make an empty file to store the final outcomes
        var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().to_array().to_numpy().squeeze()
        # var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().copy() #original, but can't s
        #Because we don't need to slice anything, we can convert it to a numpy array
        #and use Numba for faster processing
        
        obs_converted = obs_files.to_array().to_numpy().squeeze() #drop unneed dimension
        subx_converted = subx_files.to_array().to_numpy().squeeze() #drop unneed dimension
    
        
        @njit
        def skill_assessment(var_OUT, subx_converted, obs_converted):
        #Now find pearson correlation by model, lead, and lat/lon
            for model in range(var_OUT.shape[0]):
                print(f'Working on model {model+1} for pearson r correlation')
                for Y in range(var_OUT.shape[2]):
                    # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
                    for X in range(var_OUT.shape[3]):
                        for lead in np.arange(7,45,7):
                            var_OUT[model, lead,Y,X] = \
                                np.corrcoef(subx_converted[:,model, lead, Y, X],obs_converted[:,model, lead, Y, X])[0,1]
            return(var_OUT)
        
            
        #Very slow loop (I switched to numba instead)
         # #Now find pearson correlation by model, lead, and lat/lon
         # for model in range(var_OUT[f'{sub_name}'].shape[0]):
         #     print(f'Working on model {model+1} for pearson r correlation')
         #     for Y in range(var_OUT.Y.shape[0]):
         #         print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
         #         for X in range(var_OUT.X.shape[0]):
         #             for lead in np.arange(7,45,7):
         #                 var_OUT[f'{sub_name}'][model, lead,Y,X] = \
         #                     np.corrcoef(subx_files[f'{sub_name}'][:,model, lead, Y, X],obs_files[f'{sub_name}'][:,model, lead, Y, X])[0,1]
                         
        
        yearly_skill = skill_assessment(var_OUT, subx_converted, obs_converted)
        
        #Convert to an xarray object
        var_yearly = xr.Dataset(
            data_vars = dict(
                Pearson_r_coefficient = (['model','lead','Y','X'], yearly_skill[:,:,:,:]),
            ),
            coords = dict(
                X = subx_files.X.values,
                Y = subx_files.Y.values,
                lead = subx_files.lead.values,
                model = subx_files.model.values,
        
            ),
            attrs = dict(
                Description = 'Pearson r correlation over all years and lead times'),
        ) 
        
        var_yearly.to_netcdf(f'{output_nc_dir}/interannual_pearson_skill.nc')
        
    return(var_yearly)

var_yearly


def plot_interannual_pearsonR(lead_week, mod):
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.background_patch.set_facecolor('black') #changes np.nan value colors
    ax.outline_patch.set_edgecolor('black') #changes the border of the whole plot to black
    # ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k') #colors ocean blue
    
    
    ETo = var_yearly.Pearson_r_coefficient.isel(lead=lead_week).sel(model=mod)
    ETo.plot(
        transform=proj, subplot_kws={'projection':proj}, cmap='YlOrRd'
        )
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = True #this will remove the longitude lines
    # gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    plt.savefig(f'{output_image_dir}/week_lead_{lead_week//7+1}_model{mod}.tif', dpi=300)
    plt.cla()
    plt.clf()

#save plots for each model and each lead week
for mod in [1,2,3,4]:
    for lead_week in np.arange(7,49,7):
        plot_interannual_pearsonR(lead_week,mod)


'''NOW ADD MASKS TO LOOK AT REGIONS/SEASONS/AND LEADS

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

var='RZSM'
if var == 'ETo':
    obs_name = 'ETo_SubX_*.nc'
    sub_name = 'ETo_anom'
    obs_path = obs_eto_path
elif var == 'RZSM':
    obs_name = 'SM_SubX_*.nc'
    sub_name = 'RZSM_anom'
    obs_path= obs_rzsm_path

subx_files = xr.open_mfdataset(f'{var}_anomaly*.nc', concat_dim=['S'], combine='nested')
obs_files = xr.open_mfdataset(f'{obs_path}/{obs_name}', concat_dim = ['S'], combine = 'nested')

# #choose season individually 
# def summer(month):
#     return(month >=6) & (month <=8)
# def fall(month):
#     return(month >=9) & (month <=11)
# def winter(month):
#     return(month ==12) & (month ==1) & (month ==2)
# def spring(month):
#     return(month >=3) & (month <=5)

summer_subx = subx_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='JJA'))
fall_subx = subx_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='SON'))
winter_subx= subx_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='DJF'))
spring_subx= subx_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='MAM'))

summer_obs = obs_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='JJA'))
fall_obs = obs_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='SON'))
winter_obs = obs_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='DJF'))
spring_obs = obs_files[f'{sub_name}'].sel(S=(subx_files['S.season']=='MAM'))

skill_subx = [summer_subx, fall_subx, winter_subx, spring_subx]
skill_obs = [summer_obs, fall_obs, winter_obs, spring_obs]
season_name = ['Summer', 'Fall', 'Winter', 'Spring']

for season in range(len(skill_subx)):
    
    #Make an empty file to store the final outcomes
    var_OUT = skill_subx[season][0,:,:,:,:].to_dataset().to_array().to_numpy().squeeze()

# var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().copy() #original, but can't s
#Because we don't need to slice anything, we can convert it to a numpy array
#and use Numba for faster processing

    obs_converted = skill_subx[season].to_numpy().squeeze() #drop unneed dimension
    subx_converted = skill_obs[season].to_numpy().squeeze() #drop unneed dimension

       
    @njit
    def season_skill_assessment(var_OUT, subx_converted, obs_converted):
    #Now find pearson correlation by model, lead, and lat/lon
        for model in range(var_OUT.shape[0]):
            print(f'Working on model {model+1} for pearson r correlation')
            for Y in range(var_OUT.shape[2]):
                # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
                for X in range(var_OUT.shape[3]):
                    for lead in np.arange(7,45,7):
                        var_OUT[model, lead,Y,X] = \
                            np.corrcoef(subx_converted[:,model, lead, Y, X],obs_converted[:,model, lead, Y, X])[0,1]
        return(var_OUT)
    
    seasonal_skill = season_skill_assessment(var_OUT, subx_converted, obs_converted)
    
    #Convert to an xarray object
    var_yearly = xr.Dataset(
        data_vars = dict(
            Pearson_r_coefficient = (['model','lead','Y','X'], seasonal_skill[:,:,:,:]),
        ),
        coords = dict(
            X = subx_files.X.values,
            Y = subx_files.Y.values,
            lead = subx_files.lead.values,
            model = subx_files.model.values,
    
        ),
        attrs = dict(
            Description = 'Pearson r correlation for each season and lead times'),
    ) 
    
    var_yearly.to_netcdf(f'{output_nc_dir}/{season_name[season]}_pearson_skill.nc')
    
    return(var_yearly)
xr.where(a.ETo_anom == np.nan, a.ETo_anom, a.ETo_anom) = 0

'''Test with xskillscore'''
# a = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SMERGE_SM/SM_SubX_values/SM_SubX_1999-01-10.nc')
a = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/ETo_anomaly_2012-02-19.nc')
a.close()
b = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/ETo_anomaly_2012-02-19.nc')
b.close()

rzsm = xr.open_mfdataset('vas*.nc', concat_dim=['S'], combine='nested')
np.corrcoef(a.ETo_anom[0,0,7,10,10].values,a.ETo_anom[0,0,7,10,10].values)
# proj = ccrs.PlateCarree()

only_2_dim = rzsm.vas[::2,0,7,10,10].values
dim_2 = rzsm.vas[::2,0,7,11,11].values
len(only_2_dim)
np.corrcoef(only_2_dim,dim_2)
np.corrcoef(only_2_dim,dim_2)[0,1]



# create colormap 
cmap=plt.get_cmap()

# use white color to mark 'bad' values
cmap.set_bad(color='w')
# gl.xlabel_style = {'size': 15, 'color': 'gray'}
# gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
# ax.coastlines()

#example np.corrcoeff
x_simple = np.array([-2, -1, 0, 1, 2])
y_simple = np.array([4, 1, 3, 2, 0])
my_rho = np.corrcoef(x_simple, y_simple)




np.corrcoef(a.ETo_anom.isel(S=0,model=1,lead=7,Y=10,X=10).values,a.ETo_anom.isel(S=0,model=1,lead=7,Y=10,X=10).values)

ddd=np.array([1,1,1,1,2])
ccc = np.array([2,2,2,2,3])
np.corrcoef(ddd,ccc)

dataset = a
# dataset = netcdf_dataset(fname)
sst = dataset.variables['ETo_anom'][0, 0, 7,:,:]
lats = dataset.variables['Y'][:]
lons = dataset.variables['X'][:]

ax = plt.axes(projection=ccrs.PlateCarree())

plt.contourf(lons, lats, sst, 60,
             transform=ccrs.PlateCarree())

proj = ccrs.PlateCarree()


dataset.plot(transform=proj, col='Time', col_wrap=3, robust=True, subplot_kws={'projection':proj})




import cartopy.crs as ccrs
import xarray as xr




dd = xs.pearson_r(a, a, dim='S')
np.corrcoef(a.to_array(),a.to_array())

a = xr.DataArray(np.random.rand(5, 3, 3),

                 dims=['time', 'x', 'y'])

b = xr.DataArray(np.random.rand(5, 3, 3),

                 dims=['time', 'x', 'y'])

xs.pearson_r(a, b, dim='time')

xs.pearson_r(a,b)

air1d = c.isel(Y=20, X=10)
np.nanmean(air1d.to_array())
air1d.ETo_anom.plot()

plt.plot(skill['rmm1'])
plt.title('GMAO-GEOS_V2p1 RMM1 Skill')
plt.xlabel('Lead Time (Days)')
plt.ylabel('ACC')

# dask_gufunc_kwargs.setdefault("allow_rechunk", True)


HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask file
West_conus_mask = xr.open_dataset(West_conus_path)

smerge_file = xr.open_dataset(f'{dir1}/Data/SMERGE_SM/Raw_data/smerge_sm_merged_remap.nc4') #for mask

'''First, fix RZSM and ETo anomaly SubX files so that the S dimension with no data
is dropped. This will assist with attempting to do correlation using 
xr.open_mfdataset'''
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{subX_dir}/{var}*.nc4')):
        date_list.append(file[-14:-4])
    return(date_list)
        
init_date_list = return_date_list(var = 'RZSM')    


          

np.count_nonzero(open_f.where(HP_conus_mask.High_Plains > 0).to_array())
np.count_nonzero(open_f.where(HP_conus_mask.High_Plains == 0).to_array())  
    

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

#TODO: Only use the High Plains mask for the High Plains region.
#TODO: Otherwise, use the West mask for all other regions. The numbers are correct


#Test variable
var='RZSM'

os.chdir(subX_dir) #Set directory for SubX

#gridMET actual values
gmet_subX_dir = f'{dir1}/Data/gridMET/ETo_SubX_values' 

#SMERGE actual values
rzsm_subX_dir = f'{dir1}/Data/SMERGE_SM/SM_SubX_values' 
#%%
'''Now find the anomaly correlation by lead time



'''


subx_all = xr.open_mfdataset(f'{subX_dir}/{var}*delete*.nc4', concat_dim=['S'], combine='nested')
# subx_S_dim_remove =
subx_all[varname1][2001,1,:,10,20].values

obs_all = xr.open_mfdataset(f'{rzsm_subX_dir}/SM*.nc')

#lead times
start_lead = 6
lead=0
lead_weeks_as_index = np.arange(start_lead,42,7)

#Stack X and Y coordinates for easier processing
rzsm_subx_all_st = rzsm_subx_all.stack(grid = ('X','Y'))
rzsm_subx_obs_st = rzsm_obs_all.stack(grid = ('X','Y'))

#Nope
# CONUS_correlation_RZSM = xr.corr(rzsm_subx_all.RZSM_anom[0,:,start_lead+lead,:,:],rzsm_obs_all.SMERGE_SubX_value[0,:,start_lead+lead,:,:]).compute()
CONUS_correlation_RZSM = pr(rzsm_subx_all_st.RZSM_anom[0,:,start_lead+lead,:],rzsm_subx_obs_st.SMERGE_SubX_value[0,:,start_lead+lead,:])



CONUS_correlation_RZSM.unstack('grid')

#West region
week_1 = rzsm_subx_all.RZSM_anom[0,:,start_lead+lead,:,:].where(West_conus_mask == 1)

xr.corr(week_1.where(HP_conus_mask==1), week_1.where(HP_conus_mask==1)).mean().compute()

eddi_subx_all[0]

xr.corr(eddi_subx_all.where(HP_conus_mask == 1),eddi_subx_all.where(HP_conus_mask == 1)).mean().

eddi_mean = eddi_subx_all.where(HP_conus_mask == 1).mean().compute()
eddi_all = eddi_subx_all.where(HP_conus_mask == 1).compute()




#%% EDDI correlation




'''this is a test
I can only work with the exact same files because I need to re-format  SubX observations
'''

rzsm_subx_all = xr.open_mfdataset(f'{subX_dir}/{var}*.nc4', chunks={'S':-1,'model':-1,'lead':-1,'Y':-1,'X':-1}).load()

#%%


RZSM_anom_ACC = metrics.pearson_r(rzsm_subx_all, rzsm_subx_all, dim='S')
    #%%

rzsm_subx_all.sizes

rzsm_subx_all.chunks

#rechunk
rzsm_subx_all.chunk({'S':-1}).sizes
if rzsm_subx_all.chunks:
    for 
    rzsm_subx_all = rzsm_subx_all.chunk({dim:-1})

rzsm_subx_all.chunk(dict('S':-1,'model':-1,'X':-1,'Y':-1))
rzsm_obs_all = xr.open_mfdataset(f'{subX_dir}/{var}*.nc4', chunks={'S':-1})
rzsm_subx_all.RZSM_anom[0,:,:,15,15].values

#lead times
start_lead = 6
lead=0
lead_weeks_as_index = np.arange(start_lead,42,7)

#Stack X and Y coordinates for easier processing
rzsm_subx_all_st = rzsm_subx_all.stack(grid = ('X','Y'))
rzsm_subx_obs_st = rzsm_obs_all.stack(grid = ('X','Y'))

#Nope
# CONUS_correlation_RZSM = xr.corr(rzsm_subx_all.RZSM_anom[0,:,start_lead+lead,:,:],rzsm_obs_all.SMERGE_SubX_value[0,:,start_lead+lead,:,:]).compute()

# eddi_subx_all.EDDI[0,:,6,15,15].values #lead week 1

start_lead = 6
lead=0
lead_weeks_as_index = np.arange(start_lead,42,7)

week_1 = eddi_subx_all.EDDI[0,:,start_lead+lead,:,:].where(HP_conus_mask == 1)

xr.corr(week_1.where(HP_conus_mask==1), week_1.where(HP_conus_mask==1)).mean().compute()

eddi_subx_all[0]

xr.corr(eddi_subx_all.where(HP_conus_mask == 1),eddi_subx_all.where(HP_conus_mask == 1)).mean().

eddi_mean = eddi_subx_all.where(HP_conus_mask == 1).mean().compute()
eddi_all = eddi_subx_all.where(HP_conus_mask == 1).compute()















#%%
'''Since I've already created the actual values from SubX files for each index 
(EDDI, ETo, and RZSM) [from script 1d_RZSM_and_scatterplots.py], I can 
just re-factor portions of that code to calculate the pearson's correlation coefficent'''


















def flatten_RZSM_ETo_EDDI(date_list):
    '''Open RZSM anomalies and SubX file'''
    
    #Flatten 1 file to see the length
    ravel_length = len(xr.open_dataset(f'{output_dir}/ETo_SubX_{date_list[0]}.nc').to_array().values.ravel())    
    file_length = len(date_list)
    final_length = ravel_length * file_length
    
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_RZSM = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open SMERGE file and ravel, and append to empty array
        smerge_file = xr.open_dataset(f'{SM_SubX_out_dir}/SM_SubX_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        sub_file = xr.open_dataset(f'{subX_SM_out_dir}/SM_SubX_m3_m3_{date_}.nc4').SM_SubX_m3_m3_value[0,:,:,:,:].values.ravel()
        #append smerge to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_RZSM[count_value:ravel_end,0] = smerge_file[:]
        correlation_RZSM[count_value:ravel_end,1] = sub_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
 
    corr_df = pd.DataFrame(correlation_RZSM)
    corr_df.corr()
    
    '''Open gridMET and SubX file'''
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_ETo = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open gridMET file and ravel, and append to empty array
        g_file = xr.open_dataset(f'{output_dir}/ETo_SubX_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        s_file = xr.open_dataset(f'{subX_dir}/ETo_{date_}.nc4').ETo[0,:,:,:,:].values.ravel()
        #append gridMET to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_ETo[count_value:ravel_end,0] = g_file[:]
        correlation_ETo[count_value:ravel_end,1] = s_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
        
    corr_df = pd.DataFrame(correlation_ETo)
    corr_df.corr()   
    
    '''Open EDDI and SubX file'''
    #Create an array with the length final_length and with 2 columns for gridMET and SubX
    correlation_EDDI = np.empty((final_length, 2))

    count_value = 0 #Keep a counter for index
    ravel_end = ravel_length
    for date_ in date_list:
        
        #Open EDDI file and ravel, and append to empty array
        g_file = xr.open_dataset(f'{eddi_subX_dir}/EDDI_SubX_{mod}_{date_}.nc').to_array().values.ravel()
        #Open SubX ETo file. Remove 1st dimension becuase it's useless
        s_file = xr.open_dataset(f'{subX_dir}/EDDI_{date_}.nc4').EDDI[0,:,:,:,:].values.ravel()
        #append gridMET to first column (0 index)
        #append SubX to 2nd columnn (1st index)
        correlation_EDDI[count_value:ravel_end,0] = g_file[:]
        correlation_EDDI[count_value:ravel_end,1] = s_file[:]
        
        #Update counter to append properly
        count_value += ravel_length
        ravel_end += ravel_length
    
    '''Currently very low correlation'''
    corr_df = pd.DataFrame(correlation_EDDI)
    corr_df.corr()
          
    return(correlation_RZSM, correlation_ETo, correlation_EDDI)