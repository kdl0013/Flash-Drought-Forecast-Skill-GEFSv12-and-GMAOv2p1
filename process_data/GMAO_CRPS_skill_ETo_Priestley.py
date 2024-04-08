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

print(f"PYTHON: {sys.version}")  # PYTHON: 3.8.1 | packaged by conda-forge | (default, Jan 29 2020, 15:06:10) [Clang 9.0.1 ]
print(f" xarray {xr.__version__}")  # xarray 0.14.1
print(f" numpy {np.__version__}")  # numpy 1.17.3
print(f" matplotlib {mpl.__version__}")  # matplotlib 3.1.2


# dir1 = 'main_dir'


# TODO change later
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'GMAO'
name_ = 'Priestley'
# model_NAM1 = 'GMAO'
# name_='Priestley'
#%%
subX_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
os.chdir(f'{subX_dir}') #for multi ensemble mean

 #where subX model data lies 

subdir_for_mean = f'{subX_dir}/mean_for_ACC'
mask_path = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4'

#observations
obs_subx_eto_path = f'{dir1}/Data/MERRA2/ETo_SubX_values/'

#mean is needed for anomaly correlation coefficient calculation
obs_eto = xr.open_dataset(f'{dir1}/Data/MERRA2/ETo_anomaly_{name_}_MERRA.nc',chunks = {'time':1})
dates_list = list(obs_eto.time.values)
df = pd.DataFrame({'dates':dates_list})
df['dates'] = pd.to_datetime(df['dates'])
df['dates'] = df['dates'].dt.floor('d')

new_date_list = np.array(df['dates'])

obs_eto=obs_eto.assign_coords(time=new_date_list)

obs_eto_mean = xr.open_dataset(f'{dir1}/Data/MERRA2/ETo_anomaly_{name_}_MERRA_mean.nc')


output_nc_dir = f'{subX_dir}/skill_assessments/' #for skill assessment .nc4 files
output_image_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}'
output_season_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}/seasonal_skill/'

os.system(f'mkdir -p {output_season_dir}')

'''Anomaly correlation coefficient skill
Read all files in at once with xr.open_mfdataset, then calculate the skill
based only on weekly lead time. '''


#Don't have to load into memory because I'm immediately converting them into a np.array
'''ISSUE: I technically need more data from MERRA for the latest GMAO forecasts past June 20, 2022

For right now, just do through month of may'''

var='ETo'


if var == 'ETo':
    obs_name = f'ETo_SubX_anomaly_{name_}_{model_NAM1}*.nc4'
    sub_name = 'ETo_anom'

elif var == 'RZSM':
    obs_name = 'SM_SubX_anomaly_*.nc4'
    sub_name = 'RZSM_anom'

#Open files

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


subx_files = rename_subx_for_climpred(xr.open_mfdataset(f'{var}_{name_}_anomaly_{model_NAM1}_*.nc4', concat_dim=['S'], combine='nested',chunks={'S': 1, 'L': 45}).sel(S=slice('2000-01-01','2022-05-30')))
# subx_files=subx_files.chunk(dict(init=-1))
obs_eto = rename_obs_for_climpred(obs_eto)
# obs_eto = obs_eto.chunk(dict(time=-1))

fcst=subx_files.chunk({'init':-1})
verif=obs_eto.chunk({'time':-1})


#MASKS
conus_mask = xr.open_dataset(f'{mask_path}')  
HP_conus_mask = rename_obs_for_climpred(conus_mask['USDM-HP_mask'])
West_conus_mask = rename_obs_for_climpred(conus_mask['USDM-West_mask']).rename(time='init')

# #for GMAO and RSMAS
# #Add a mask
# for i in range(fcst.init.shape[0]):
#     single_day = fcst.ETo_anom[i:i+1,:,:,:,:]
#     single_day.where(West_conus_mask[0,:,:]==1)
#     out_ = xr.concat(fcst.ETo_anom[i:i+1,:,:,:,:])
#     print(i)
# fcst.where(West_conus_mask==1)

#%%
# #ACC skill
# skill = hindcast.verify(
#     metric="acc", comparison="e2o", dim="init", alignment="maximize"
# )
# skill
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


def select_region(cluster_num):
    
    output_dictionary = {}
   
    #restrict season
    for season in ['JJA','SON','DJF','MAM']:

        skill_season=skillds.sel(init=(skillds['init.season']==f'{season}'))
    
        # skillds.ETo_anom.shape
        if cluster_num == 1:
            all_crpss_skill = skill_season.ETo_anom.where(West_conus_mask[0,:,:]== 1)
            # bn.nanmean(all_crpss_skill)
        else:
            all_crpss_skill = skill_season.ETo_anom.where(West_conus_mask[0,:,:] == cluster_num)
        
        #Keep only weekly leads
        skill_lead=all_crpss_skill[::7,:,:,:]
        skill_lead=skill_lead[1:,:,:,:].to_dataset()
        
        #lead_vals = 
        for lead in range(skill_lead.lead.shape[0]):
            
            output_dictionary[f'Lead_{lead+1}_{season}']= skill_lead.ETo_anom[lead,:,:,:].mean().values
        
    return(output_dictionary)
#%%
#add all datasets together
hindcast = HindcastEnsemble(fcst).add_observations(verif)
skillds = hindcast.verify(metric="crpss", comparison="m2o", dim="member", alignment="maximize").persist()
#%%

CRPSS_skill = {}
for clus_num in np.arange(1,7):
    CRPSS_skill[f'Cluster {clus_num}'] = select_region(cluster_num=clus_num)
    #print(all_cluster_acc_ETo[f'Cluster {clus_num}'])

#%%
def setup_plot_all_leads(CRPSS_skill,cluster_num):
    #The goal of this function is to take our outputs from anomaly correlation 
    #coefficient, convert them to a dataframe by model, lead, season,
    #to set up for plot
    col_names = ['Lead','Season','CRPSS']
    
    df_OUT = pd.DataFrame(columns=col_names)
    
    var_cluster = CRPSS_skill[f'Cluster {cluster_num}']
    
    for idx,k in enumerate(var_cluster.keys()):
        # if idx==1:
        #     break
        #split to get model, lead, and season
        split_ = k.split('_')
        
        df_OUT=df_OUT.append({'Lead':split_[1], \
                             'Season': split_[-1],'CRPSS':var_cluster[k]},ignore_index=True)

    return(df_OUT)

all_vals_setup = {}
for clus_num in np.arange(1,7):
    all_vals_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(CRPSS_skill=CRPSS_skill,cluster_num=clus_num)

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
            var_OUT[0,:] = subset_by_cluster[subset_by_cluster['Season']=='MAM']['CRPSS']
            var_OUT[1,:] = subset_by_cluster[subset_by_cluster['Season']=='JJA']['CRPSS']
            var_OUT[2,:] = subset_by_cluster[subset_by_cluster['Season']=='SON']['CRPSS']
            var_OUT[3,:] = subset_by_cluster[subset_by_cluster['Season']=='DJF']['CRPSS']
            
            #Set up for p_color mesh
            x=np.arange(1.5,var_OUT.shape[1]+2)
            x = x-1
            y=np.arange(1.5,var_OUT.shape[0]+2)
            Z=var_OUT[:,:]
            
            # #now only grab the weekly leads
            # var_LEAD = var_OUT[:,::7]
            # var_LEAD = var_LEAD[:,1:]
            
            Index = ['Spring', 'Summer', 'Fall', 'Winter']
            if model_NAM1 == 'ESRL' or model_NAM1 == 'EMC':
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
 
    cmap = mpl.colors.ListedColormap(plt.cm.RdPu(np.linspace(min_all, max_all, 7)))
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
                s.set_title(f'RSMAS CCSM4 \n ETo {name_} Continous Ranked Probability Skill Score',fontsize=25)
            elif model_NAM1 == 'GMAO':
                s.set_title(f'GMAO GEOS-5 \n ETo {name_} Continous Ranked Probability Skill Score',fontsize=25)
            elif model_NAM1 == 'ESRL':
                s.set_title(f'ESRL FIMr1p1 \n ETo {name_} Continous Ranked Probability Skill Score',fontsize=25)
            elif model_NAM1 == 'EMC':
                s.set_title(f'EMC GEFSv12 \n ETo {name_} Continous Ranked Probability Skill Score',fontsize=25)

    s.set_xlabel('Week Lead',fontsize=25)
    plt.savefig(f'{output_season_dir}/all_season_CRPSS_skill_{var}_{name_}.tif',dpi=300)
        
    return(0)

make_save_plots_all_models_seasons_leads(acc_values_to_plot=acc_values_to_plot,var=var,min_all=min_all,max_all=max_all)


