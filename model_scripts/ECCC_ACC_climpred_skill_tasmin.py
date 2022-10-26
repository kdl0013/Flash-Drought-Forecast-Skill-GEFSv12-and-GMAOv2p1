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

NOW ADD MASKS TO LOOK AT REGIONS/SEASONS/AND LEADS

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

TODO: Use only the West mask for the West region. Use HP_conus_mask for all other regions.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
import numpy as np
import matplotlib as mpl
import os
import pandas as pd
from glob import glob
import sys
from numpy import inf
from multiprocessing import Pool
from os.path import exists
from datetime import datetime
from climpred import HindcastEnsemble
from tqdm.auto import tqdm
import bottleneck as bn
from climpred import metrics as cmet
import dask.dataframe as dd
import datetime as dt
import dask
from dask import array
import pickle
from os.path import exists
print(f"PYTHON: {sys.version}")  # PYTHON: 3.8.1 | packaged by conda-forge | (default, Jan 29 2020, 15:06:10) [Clang 9.0.1 ]
print(f" xarray {xr.__version__}")  # xarray 0.14.1
print(f" numpy {np.__version__}")  # numpy 1.17.3
print(f" matplotlib {mpl.__version__}")  # matplotlib 3.1.2

dask.config.set(**{'array.slicing.split_large_chunks': False})
# dir1 = 'main_dir'


# TODO change later
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'ECCC'
var = 'tasmin'
# model_NAM1 = 'ESRL'
# name_='Priestley'
subX_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
os.chdir(f'{subX_dir}') #for multi ensemble mean

output_nc_dir = f'{subX_dir}/skill_assessments/' #for skill assessment .nc4 files

file_exists = exists(f'{output_nc_dir}/{model_NAM1}_{var}_anomaly_climpred_ACC.pkl')

if file_exists:
    print(f'Completed {model_NAM1} {var}.')
    sys.exit(0)

#%%


 #where subX model data lies 
subdir_for_mean = f'{subX_dir}/mean_for_ACC'
mask_path = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4'

#observations
obs_subx_eto_path = f'{dir1}/Data/MERRA2/{var}_SubX_values/'

#return name of xarray variable
def name(file):
    return(list(file.keys())[0])

#mean is needed for anomaly correlation coefficient calculation

if var == 'srad':
    obs_eto = xr.open_dataset(f'{dir1}/Data/MERRA2/radiation_anomaly_MERRA.nc',chunks = {'time':1})
    obs_eto_mean = xr.open_dataset(f'{dir1}/Data/MERRA2/radiation_anomaly_MERRA_mean.nc')
    
    #check values
    # print(obs_eto[name(obs_eto)][0,10,10].values)
    # print(obs_eto_mean[name(obs_eto_mean)][0,10,10].values)
    # print(np.nanmax(obs_eto[name(obs_eto)][:,10,10].values))
    
else:
    obs_eto = xr.open_dataset(f'{dir1}/Data/MERRA2/{var}_anomaly_MERRA.nc',chunks = {'time':1})
    obs_eto_mean = xr.open_dataset(f'{dir1}/Data/MERRA2/{var}_anomaly_MERRA_mean.nc')
    # print(obs_eto[name(obs_eto)][0,10,10].values)
    # print(obs_eto_mean[name(obs_eto_mean)][0,10,10].values)
    # print(np.nanmax(obs_eto[name(obs_eto)][:,10,10].values))
    
dates_list = list(obs_eto.time.values)
df = pd.DataFrame({'dates':dates_list})
df['dates'] = pd.to_datetime(df['dates'])
df['dates'] = df['dates'].dt.floor('d')

new_date_list = np.array(df['dates'])

obs_eto=obs_eto.assign_coords(time=new_date_list)

output_image_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}'
output_season_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}/seasonal_skill/'

os.system(f'mkdir -p {output_season_dir}')

'''Anomaly correlation coefficient skill
Read all files in at once with xr.open_mfdataset, then calculate the skill
based only on weekly lead time. '''

#Don't have to load into memory because I'm immediately converting them into a np.array
'''ISSUE: I technically need more data from MERRA for the latest GMAO forecasts past June 20, 2022

For right now, just do through month of may'''

metric='ACC'
obs_name = f'{var}_SubX_anomaly_{model_NAM1}*.nc4'


#Open files
#%%

def rename_subx_for_climpred(file):
    #https://climpred.readthedocs.io/en/stable/examples/subseasonal/daily-subx-example.html
    file = file.rename(S='init')
    file = file.rename(L='lead')
    file["lead"].attrs = {"units": "days"}
    file = file.rename(M='member')
    file = file.rename(X='lon')
    file = file.rename(Y='lat')
    file = file.assign_attrs(lead='days')
    return(file)

def rename_obs_for_climpred(file):
    file = file.rename(X='lon')
    file = file.rename(Y='lat')
    return(file)

if model_NAM1 == 'ESRL':
    # subx_files = subx_files.sel(S=obs_files.S.values)
    #for some reason 5 dates aren't correct, fix them (doesn't work)
    # for i in sorted(glob(f'{var}_{name_}_anomaly_{model_NAM1}_*.nc4')):
    #     open_f = xr.open_dataset(i)
    #     open_f.close()
    #     open_f=open_f.assign_coords(S=np.atleast_1d(pd.to_datetime(i.split('_')[-1].split('.')[0])))
    #     open_f.to_netcdf(f'{var}_{name_}_anomaly_{model_NAM1}_*.nc4')
        
    subx_files = rename_subx_for_climpred(xr.open_mfdataset(f'{var}_anomaly_{model_NAM1}_*.nc4', concat_dim=['S'], combine='nested',chunks={'S': 1, 'L': 32}).sel(S=slice('2000-01-01','2022-05-30')))
    subx_files=subx_files.sel(init=~subx_files.get_index('init').duplicated())
    '''ISSUES with only this dataset. Have no idea why, it creates a duplicate. And 5 files (which I have checked have the right names,
    those five files will appear with dates of 0'''
    #Just create a new list of datetime dates
    start_ = dt.date(2000,1,1)
    end_ = dt.date(2022,5,25)
    #dates
    dates = [start_ + dt.timedelta(days=d) for d in range(0, end_.toordinal() - start_.toordinal() + 1)]
    ###only get wednesday initialized days
    dates = [pd.to_datetime(i) for i in dates if i.weekday() == 2]
    # dates[0]
    len(dates)
    
    dates_ = [i for i in dates if i in subx_files.init.values ]
    subx_files=subx_files.assign_coords(init=dates_)

else:
    subx_files = rename_subx_for_climpred(xr.open_mfdataset(f'{var}_anomaly_{model_NAM1}_*.nc4', concat_dim=['S'], combine='nested',chunks={'S': 1}).sel(S=slice('2000-01-01','2022-05-30')))
    
if var == 'tas':
    #for some reason there is an issue with the dates, maybe try tasmax which already works
    test_tasmax = rename_subx_for_climpred(xr.open_mfdataset(f'tasmax_anomaly_{model_NAM1}_*.nc4', concat_dim=['S'], combine='nested',chunks={'S': 1}).sel(S=slice('2000-01-01','2022-05-30')))
    subx_files['init'] = test_tasmax.init.values
# subx_files=subx_files.chunk(dict(init=-1))
obs_eto = rename_obs_for_climpred(obs_eto)
# obs_eto = obs_eto.chunk(dict(time=-1))

#Just rename to help out climpred
fcst=subx_files.chunk({'init':-1}).rename({name(subx_files) : f'{var}'})
# check_vals = (fcst[name(fcst)][:,0,:,10,10].values)
# np.nanmax(check_vals)

verif=obs_eto.chunk({'time':-1}).rename({name(obs_eto) : f'{var}'})


#MASKS
conus_mask = xr.open_dataset(f'{mask_path}')  
HP_conus_mask = rename_obs_for_climpred(conus_mask['USDM-HP_mask']).rename(time='init')
West_conus_mask = rename_obs_for_climpred(conus_mask['USDM-West_mask']).rename(time='init')

#%%
#TODO: Add each region and save
def return_cluster_name(cluster_num):
    if cluster_num == 1:
        out_name= 'West' 
    elif cluster_num == 2:
        out_name= 'Midwest'
    elif cluster_num == 3:
        out_name= 'HighPlains'        
    elif cluster_num == 4:
        out_name= 'South'      
    elif cluster_num == 5:
        out_name= 'Southeast'
    elif cluster_num == 6:
        out_name= 'Northeast'  
    return(out_name)

#%%
def select_region(cluster_num):
    output_dictionary = {}
    #restrict season
    for season in ['JJA','SON','DJF','MAM']:
     
        # restrict to only a season first, then compute the ACC score
        #for some reason climpred always removes the 
        if cluster_num == 1:
            hindcast = HindcastEnsemble(fcst[name(fcst)].where(West_conus_mask[0,:,:]== 1).sel(init=(fcst['init.season']==f'{season}'))).add_observations(verif[name(verif)].where(West_conus_mask[0,:,:]== 1))

            # bn.nanmean(all_crpss_skill)
        else:
            #Keep all seasons
            hindcast = HindcastEnsemble(fcst[name(fcst)].where(HP_conus_mask[0,:,:]== cluster_num).sel(init=(fcst['init.season']==f'{season}'))).add_observations(verif[name(verif)].where(West_conus_mask[0,:,:]== cluster_num))
        
        skillACC = hindcast.verify(metric=f"{metric.lower()}", comparison="e2o", dim="init", alignment="maximize")
        # bn.nanmean(skillACC.ETo_anom)
        skill_lead=skillACC[name(skillACC)][::7,:,:]
        skill_lead=skill_lead[1:,:,:] #don't want the 0th day (initialization day)
        # skill_lead=skill_lead[1:,:,:].compute().to_dataset()
        #lead_vals = 
        for lead in range(skill_lead.lead.shape[0]):
            output_dictionary[f'Lead_{lead+1}_{season}']= skill_lead[lead,:,:].mean().values

        
    return(output_dictionary)




ACC_skill  = {}
for clus_num in np.arange(1,7):
    print(f'Working on Cluster {clus_num} for variable {var} model {model_NAM1}.')
    ACC_skill[f'Cluster {clus_num}'] = select_region(cluster_num=clus_num)

#%%
def setup_plot_all_leads(ACC_skill,cluster_num):
    #The goal of this function is to take our outputs from anomaly correlation 
    #coefficient, convert them to a dataframe by model, lead, season,
    #to set up for plot
    col_names = ['Lead','Season',f'{metric}']
    
    df_OUT = pd.DataFrame(columns=col_names)
    
    var_cluster = ACC_skill[f'Cluster {cluster_num}']
    
    for idx,k in enumerate(var_cluster.keys()):
        # if idx==1:
        #     break
        #split to get model, lead, and season
        split_ = k.split('_')
        
        df_OUT=df_OUT.append({'Lead':split_[1], \
                             'Season': split_[-1],f'{metric}':var_cluster[k]},ignore_index=True)

    return(df_OUT)

all_vals_setup = {}
for clus_num in np.arange(1,7):
    all_vals_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(ACC_skill=ACC_skill,cluster_num=clus_num)

#%%
def plot_lead_week_season_model(all_vals_setup):
    
    #for each cluster, get all the data prepared for heatmaps
    out_clusters = {}  
    for i in all_vals_setup.keys():
        # all_vals_setup[i]

        cluster_num = int(i.split()[-1])
          
        out_name = return_cluster_name(cluster_num)                
    
        def return_data_for_plot(subset_by_cluster):
            #convert to np array for pcolor mesh
            #season order: Spring, Summer, Fall, Winter
            
            length_of_leads = (all_vals_setup[list(all_vals_setup.keys())[0]]['Season'] == 'DJF').sum()
            var_OUT = np.zeros(shape = (4,length_of_leads))
            #Place names in this order for better visual
            var_OUT[0,:] = subset_by_cluster[subset_by_cluster['Season']=='MAM'][f'{metric}']
            var_OUT[1,:] = subset_by_cluster[subset_by_cluster['Season']=='JJA'][f'{metric}']
            var_OUT[2,:] = subset_by_cluster[subset_by_cluster['Season']=='SON'][f'{metric}']
            var_OUT[3,:] = subset_by_cluster[subset_by_cluster['Season']=='DJF'][f'{metric}']
            
            #Set up for p_color mesh
            x=np.arange(1.5,var_OUT.shape[1]+2)
            x = x-1
            y=np.arange(1.5,var_OUT.shape[0]+2)
            Z=var_OUT[:,:]
            
            # #now only grab the weekly leads
            # var_LEAD = var_OUT[:,::7]
            # var_LEAD = var_LEAD[:,1:]
            
            Index = ['Spring', 'Summer', 'Fall', 'Winter']
            if model_NAM1 == 'ESRL' or model_NAM1 == 'EMC' or model_NAM1 == 'ECCC':
                Cols = ['1', '2', '3', '4']
            else:
                Cols = ['1', '2', '3', '4', '5', '6']
            var_OUT = pd.DataFrame(var_OUT,index=Index, columns=Cols)
            
            return(var_OUT,x,y,Z)
        
        #p0-p3 contains the necessary data in a dataframe
        #p0x...p0y etc just contain data that was used for old plotting, not needed anymore
        p0,p0x,p0y,p0z = return_data_for_plot(subset_by_cluster=all_vals_setup[i])

        # multi_mean = (p0+p1+p2+p3)/4
        
        # out_arrays = [p0,p1,p2,p3]
        # #add all data to a dictionary

        out_clusters[f'{out_name}'] = p0
          
        #add multi_mean to dictionary
        # out_clusters[f'{out_name}_mean'] = multi_mean
            
    return(out_clusters)

acc_values_to_plot= plot_lead_week_season_model(all_vals_setup=all_vals_setup)

def get_max_and_min_all_clusters(acc_values_to_plot):
    max_ = 0
    min_ = 1
    
    for i in acc_values_to_plot.keys():
        if np.max(acc_values_to_plot[i].values) > max_:
            max_ = np.max(acc_values_to_plot[i].values)
        if np.min(acc_values_to_plot[i].values) < min_:
            min_ = np.min(acc_values_to_plot[i].values)
    return(max_, min_)

max_all, min_all = get_max_and_min_all_clusters(acc_values_to_plot)     
#%%
def make_save_plots_all_models_seasons_leads(acc_values_to_plot,var,min_all,max_all):
    parameters = {'axes.labelsize': 16,
              'axes.titlesize': 20,
              'xtick.labelsize' : 15,
              'figure.titlesize': 24}
    plt.rcParams.update(parameters)
 
    # create figure
 
    cmap = mpl.colors.ListedColormap(plt.cm.Reds(np.linspace(min_all, max_all, 7)))
    # cmap.set_under((.8, .8, .8, 1.0))
 
    fig, ax = plt.subplots(6, 1, figsize=(29, 10), dpi=300)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    
    region_name = ''
    region_index= 0
    count_index = 0
    
    for idx,region_model in enumerate(sorted(acc_values_to_plot)):
        region_current = region_model.split('_')[0] #get region name

        if idx == len(sorted(acc_values_to_plot))-1:
            s=sns.heatmap(ax=ax[idx],data = acc_values_to_plot[region_current],
                        cmap=cmap, vmin=min_all, vmax=max_all, annot=True, 
                        fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False,cbar_ax= cbar_ax,
                        xticklabels=True,cbar=True)
        else:
            s=sns.heatmap(ax=ax[idx],data = acc_values_to_plot[region_current],
                        cmap=cmap, vmin=min_all, vmax=max_all, annot=True, 
                        fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False,cbar_ax= cbar_ax,
                        xticklabels=False,cbar=True)
        

        
        ax[idx].set_ylabel(region_current)
        if idx==0:
            if model_NAM1 == 'RSMAS':
                s.set_title(f'RSMAS CCSM4 \n {var} anomaly {metric}',fontsize=25)
            elif model_NAM1 == 'GMAO':
                s.set_title(f'GMAO GEOS-5 \n {var} anomaly {metric}',fontsize=25)
            elif model_NAM1 == 'ESRL':
                s.set_title(f'ESRL FIMr1p1 \n {var} anomaly {metric}',fontsize=25)
            elif model_NAM1 == 'EMC':
                s.set_title(f'EMC GEFSv12 \n {var} anomaly {metric}',fontsize=25)
            elif model_NAM1 == 'ECCC':
                s.set_title(f'ESRL FIMr1p1 \n {var} anomaly {metric}',fontsize=25)
            elif model_NAM1 == 'NRL':
                s.set_title(f'EMC GEFSv12 \n {var} anomaly {metric}',fontsize=25)

    s.set_xlabel('Week Lead',fontsize=25)
    plt.savefig(f'{output_season_dir}/all_season_{metric}_skill_{var}_anomaly_climpred.tif',dpi=300)
        
    return(0)

make_save_plots_all_models_seasons_leads(acc_values_to_plot=acc_values_to_plot,var=var,min_all=min_all,max_all=max_all)

#TODO save data in a file to not have to re-run anything
os.system(f'mkdir -p {output_nc_dir}')

acc_values_to_plot
f = open(f'{output_nc_dir}/{model_NAM1}_{var}_anomaly_climpred_ACC.pkl','wb')
pickle.dump(acc_values_to_plot,f)
f.close()

