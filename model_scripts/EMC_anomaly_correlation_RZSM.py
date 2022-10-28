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
import datetime

print(f"PYTHON: {sys.version}")  # PYTHON: 3.8.1 | packaged by conda-forge | (default, Jan 29 2020, 15:06:10) [Clang 9.0.1 ]
print(f" xarray {xr.__version__}")  # xarray 0.14.1
print(f" numpy {np.__version__}")  # numpy 1.17.3
print(f" matplotlib {mpl.__version__}")  # matplotlib 3.1.2


# dir1 = 'main_dir'


# TODO change later
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'EMC'

# model_NAM1 = 'RSMAS'
#%%
subX_dir = f'{dir1}/Data/SubX/{model_NAM1}/anomaly'
os.chdir(f'{subX_dir}/MEM') #for multi ensemble mean

 #where subX model data lies 

subdir_for_mean = f'{subX_dir}/mean_for_ACC'
mask_path = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4'

#observations
obs_subx_eto_path = f'{dir1}/Data/MERRA2/RZSM_SubX_values/'

#mean is needed for anomaly correlation coefficient calculation
obs_eto_mean = xr.open_dataset(f'{dir1}/Data/MERRA2/RZSM_anomaly_MERRA.nc')

output_nc_dir = f'{subX_dir}/skill_assessments/' #for skill assessment .nc4 files
output_image_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}'
output_season_dir = f'{dir1}/Outputs/anomaly_correlation/{model_NAM1}/seasonal_skill/'

os.system(f'mkdir -p {output_season_dir}')

'''Anomaly correlation coefficient skill
Read all files in at once with xr.open_mfdataset, then calculate the skill
based only on weekly lead time. '''


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


#Don't have to load into memory because I'm immediately converting them into a np.array
'''ISSUE: I technically need more data from MERRA for the latest GMAO forecasts past June 20, 2022

For right now, just do through month of may'''

var='RZSM'

conus_mask = xr.open_dataset(f'{mask_path}')  
HP_conus_mask = conus_mask['USDM-HP_mask']
West_conus_mask = conus_mask['USDM-West_mask']


obs_name = 'RZSM_SubX_anomaly_*.nc4'
sub_name = 'RZSM_anom'
sub_name_MEM = 'RZSM_anom_MEM'

#Open files
subx_files = xr.open_mfdataset(f'{var}_anomaly_MEM_{model_NAM1}_*.nc4', concat_dim=['S'], combine='nested').sel(S=slice('2000-01-01','2022-05-30'))   
obs_files = xr.open_mfdataset(f'{obs_subx_eto_path}/{obs_name}', concat_dim = ['S'], combine = 'nested').sel(S=slice('2000-01-01','2022-05-30'))    

if model_NAM1 == 'RSMAS':
    all_dates = [i for i in obs_files.S.values]
    subx_files=subx_files.sel(S=~subx_files.get_index('S').duplicated()) #remove duplicate that was causing an issue
    
else:
    subx_files = subx_files.sel(S=obs_files.S.values)
#%%
#variables for function (test)

# cluster_num = 1
# var='RZSM'
def all_season_mod_skill(cluster_num,obs_files,subx_files):
    '''Split by season, model, lead, then calculate the anomaly correlation 
    coefficient'''
    print(f'Working on cluster {cluster_num} out of 6 for RZSM.')

    # subx_files.ETo_anom[0,3,7,10,10].values
    # obs_files.ETo_anom[0,3,7,10,10].values
    ''''Reassign coordinates because of a random issue where if you don't fix it
    there is one day missing when running the next block of code between subx and obs for each season'''
    obs_files=obs_files.assign_coords(S=subx_files.S.values)

    '''PICK SEASON and REGION'''
    #only choose West conus mask file for west cluster
    if cluster_num == 1:
        summer_subx = subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='JJA')).where(West_conus_mask == cluster_num)
        
        fall_subx = subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='SON')).where(West_conus_mask == cluster_num)
        winter_subx= subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='DJF')).where(West_conus_mask == cluster_num)
        spring_subx= subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='MAM')).where(West_conus_mask == cluster_num)
        
        summer_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='JJA')).where(West_conus_mask == cluster_num)
        fall_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='SON')).where(West_conus_mask == cluster_num)
        winter_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='DJF')).where(West_conus_mask == cluster_num)
        spring_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='MAM')).where(West_conus_mask == cluster_num)
        
    elif cluster_num >=2 and cluster_num <=6:
        summer_subx = subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='JJA')).where(HP_conus_mask == cluster_num)
        fall_subx = subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='SON')).where(HP_conus_mask == cluster_num)
        winter_subx= subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='DJF')).where(HP_conus_mask == cluster_num)
        spring_subx= subx_files[f'{sub_name_MEM}'].sel(S=(subx_files['S.season']=='MAM')).where(HP_conus_mask == cluster_num)
        
        summer_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='JJA')).where(HP_conus_mask == cluster_num)
        fall_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='SON')).where(HP_conus_mask == cluster_num)
        winter_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='DJF')).where(HP_conus_mask == cluster_num)
        spring_obs = obs_files[f'{sub_name}'].sel(S=(obs_files['S.season']=='MAM')).where(HP_conus_mask == cluster_num)
    
    #make lists to iterate through them
    skill_subx = [summer_subx, fall_subx, winter_subx, spring_subx]
    skill_obs = [summer_obs, fall_obs, winter_obs, spring_obs]
    season_name = ['Summer', 'Fall', 'Winter', 'Spring']
    
    #save outputs for each model, lead, season
    output_dictionary = {}
    for season in range(len(skill_subx)):
        subx_converted = (skill_subx[season].to_numpy().squeeze()) #drop unneed dimension
        # subx_converted[:,:,10,10]
        subx_converted=subx_converted[:,::7,:,:]
        subx_converted=subx_converted[:,1:,:,:]
        # subx_converted[:,:,10,10]
        obs_converted = (skill_obs[season].to_numpy().squeeze())  #drop unneed dimension
        obs_converted=obs_converted[:,::7,:,:]
        obs_converted=obs_converted[:,1:,:,:]
        obs_converted[:,:,10,10]
        # obs_converted = obs_converted[:,0,:,:,:] #only need 1 set of observations because they are all the same
        
        #Make an empty file to store the fcluster_numinal outcomes
        if model_NAM1 == "GMAO" or model_NAM1 == "RSMAS":
            add_ = 4
        else:
            add_ = 1
        var_OUT = np.zeros(shape=(subx_converted.shape[0],(subx_converted.shape[1]+add_),subx_converted.shape[2],subx_converted.shape[3]))
        var_OUT[:,:,:,:] = np.nan
        #Only want to save dim (lead x Y x X)
        var_OUT = var_OUT[0,:,:,:]

        def season_anomaly_correlation_coefficient(var_OUT, subx_converted, obs_converted):

            '''I put this function into the loop (tried numba, didn't work well) 
            Source ACC:
            https://metclim.ucd.ie/wp-content/uploads/2017/07/DeterministicSkillScore.pdf
            def ACC(FC_anom,OBS_anom):
                top = np.nanmean(FC_anom*OBS_anom) #all forecast anomalies * all observation anomalies                    
                bottom = np.sqrt(np.nanmean(FC_anom**2)*np.nanmean(OBS_anom**2)) #variance of forecast anomalies * variance of observation anomalies
                ACC = top/bottom
                return (ACC)
            '''
            
            # #Now find pearson correlation by model, lead, and lat/lon
            # for model in range(var_OUT.shape[0]):
            # # for model in range(3,4): for testing
            #     # print(f'Working on model {model+1} for pearson r correlation')
            for Y in range(var_OUT.shape[1]):
                # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
                for X in range(var_OUT.shape[2]):
                    #anomalies are only present in weekly leads of 7
                    for lead in range(obs_converted.shape[1]):
                        
                        '''There is a Zero division error that occurs, to fix this (because numba doesn't like it)
                        just check and see if the two files have all 0s or np.nans'''

                        #ACC from function
                        obs = obs_converted[:, lead, Y, X]
                        subx = subx_converted[:, lead, Y, X]
                        
                        top = np.nanmean(subx*obs)
             
                        bottom = np.sqrt(np.nanmean(subx**2)*np.nanmean(obs**2))
                        ACC = top/bottom
                        
                        var_OUT[lead,Y,X] = ACC

            return(var_OUT)
        
        seasonal_skill = season_anomaly_correlation_coefficient(var_OUT, subx_converted, obs_converted)
        # seasonal_skill=var_OUT
        
        # #Now add back to a dictionary for each model/lead week/weason
        # for lead_week in range(obs_converted.shape[1]):
        #     seasonal_mod_skill = np.nanmean(seasonal_skill[lead_week,:,:])
        #     output_dictionary[f'Lead{lead_week+1}_{season_name[season]}']= seasonal_mod_skill
    
    
        #TODO: Now get the averages of weekly leads
        def weekly_averages(var_OUT):
            
            '''I put this function into the loop (tried numba, didn't work well) 
            Source ACC:
            https://metclim.ucd.ie/wp-content/uploads/2017/07/DeterministicSkillScore.pdf
            def ACC(FC_anom,OBS_anom):
                top = np.nanmean(FC_anom*OBS_anom) #all forecast anomalies * all observation anomalies                    
                bottom = np.sqrt(np.nanmean(FC_anom**2)*np.nanmean(OBS_anom**2)) #variance of forecast anomalies * variance of observation anomalies
                ACC = top/bottom
                return (ACC)
            '''


            # #Now find pearson correlation by model, lead, and lat/lon
            # for model in range(var_OUT.shape[0]):
            # # for model in range(3,4): for testing
            #     # print(f'Working on model {model+1} for pearson r correlation')
            for Y in range(var_OUT.shape[1]):
                # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
                for X in range(var_OUT.shape[2]):
                    
                    #week 3-4
                    if model_NAM1 == 'GMAO' or model_NAM1 == 'RSMAS':

                        wk34=obs_converted[:, 2:4, Y, X].mean(axis=1)
                        wk345=obs_converted[:, 2:5, Y, X].mean(axis=1)
                        wk3456=obs_converted[:, 2:, Y, X].mean(axis=1)
                        wk456=obs_converted[:, 3:, Y, X].mean(axis=1)
                        
                        bwk34=subx_converted[:, 2:4, Y, X].mean(axis=1)
                        bwk345=subx_converted[:, 2:5, Y, X].mean(axis=1)
                        bwk3456=subx_converted[:, 2:, Y, X].mean(axis=1)
                        bwk456=subx_converted[:, 3:, Y, X].mean(axis=1)
                        
                        
                        top34 = np.nanmean(wk34*bwk34)
                        top345 = np.nanmean(wk345*bwk345)
                        top3456 = np.nanmean(wk3456*bwk3456)
                        top456 = np.nanmean(wk456*bwk456)
                    
                        bottom34 = np.sqrt(np.nanmean(wk34**2)*np.nanmean(bwk34**2))
                        bottom345 = np.sqrt(np.nanmean(wk345**2)*np.nanmean(bwk345**2))
                        bottom3456 = np.sqrt(np.nanmean(wk3456**2)*np.nanmean(bwk3456**2))
                        bottom456 = np.sqrt(np.nanmean(wk456**2)*np.nanmean(bwk456**2))
                        
                        
                        ACC34 = top34/bottom34
                        ACC345 = top345/bottom345
                        ACC3456 = top3456/bottom3456
                        ACC456 = top456/bottom456
                        
                        var_OUT[-4,Y,X] = ACC34
                        var_OUT[-3,Y,X] = ACC345
                        var_OUT[-2,Y,X] = ACC3456
                        var_OUT[-1,Y,X] = ACC456

                    else:
                        obs_converted.shape
                        wk34=obs_converted[:, 2:4, Y, X].mean(axis=1)
                        bwk34=subx_converted[:, 2:4, Y, X].mean(axis=1)
                        top34 = np.nanmean(wk34*bwk34)
                        bottom34 = np.sqrt(np.nanmean(wk34**2)*np.nanmean(bwk34**2))
                        ACC34 = top34/bottom34
                        var_OUT[-1,Y,X] = ACC34

            return(var_OUT)

        seasonal_skill = weekly_averages(seasonal_skill)
        # seasonal_skill=var_OUT
        
        #Now add back to a dictionary for each model/lead week/weason
        for lead_week in range(seasonal_skill.shape[0]):
            seasonal_mod_skill = np.nanmean(seasonal_skill[lead_week,:,:])
            if model_NAM1=='GMAO' or model_NAM1 == 'RSMAS':
                if lead_week == 6:
                    lead_week_='3.4'
                elif lead_week == 7:
                    lead_week_='3.5'
                elif lead_week == 8:
                    lead_week_='3.6'
                elif lead_week == 9:
                    lead_week_='4.6'
                else:
                    lead_week_ = str(lead_week+ 1) 
                # print(lead_week_)
                output_dictionary[f'Lead_{lead_week_}_{season_name[season]}']= seasonal_mod_skill
            else :
                if lead_week == 4:
                    lead_week_='3.4'
                else:
                    lead_week_ = str(lead_week+ 1) 
                # print(lead_week_)
                output_dictionary[f'Lead_{lead_week_}_{season_name[season]}']= seasonal_mod_skill
    

    # all_cluster_acc_ETo[f'Cluster {clus_num}'] = output_dictionary
    
    #just appending to a dictionary
    
    return(output_dictionary)

#%%
'''mulitprocessing doesn't seem to work correctly in this situation'''
# if __name__ == '__main__':
#     p = Pool(6)
#     global all_cluster_acc_ETo
#     all_cluster_acc_ETo = {}
#     p.map(all_season_mod_skill,np.array(np.arange(1,7)))



all_cluster_acc_ETo = {}
for clus_num in np.arange(1,7):
    all_cluster_acc_ETo[f'Cluster {clus_num}'] = all_season_mod_skill(cluster_num=clus_num,obs_files=obs_files,subx_files=subx_files)
    #print(all_cluster_acc_ETo[f'Cluster {clus_num}'])

#%% Plot anomaly values in pcolormesh

#test
# var_cluster=output_dictionary
# r1_rzsm=output_dictionary

def setup_plot_all_leads(all_cluster_acc,cluster_num):
    #The goal of this function is to take our outputs from anomaly correlation 
    #coefficient, convert them to a dataframe by model, lead, season,
    #to set up for plot
    col_names = ['Lead','Season','ACC']
    
    df_OUT = pd.DataFrame(columns=col_names)
    
    var_cluster = all_cluster_acc[f'Cluster {cluster_num}']
    
    for idx,k in enumerate(var_cluster.keys()):
        # if idx==1:
        #     break
        #split to get model, lead, and season
        split_ = k.split('_')
        
        df_OUT=df_OUT.append({'Lead':split_[1], \
                             'Season': split_[-1],'ACC':var_cluster[k]},ignore_index=True)

    return(df_OUT)

# def setup_plot_mean_of_multi_weeks(all_cluster_acc,cluster_num):
#     #The goal of this function is to take our outputs from anomaly correlation 
#     #coefficient, convert them to a dataframe by model, lead, season,
#     #to set up for plot
#     col_names = ['Lead','Season','ACC']
    
#     df_OUT = pd.DataFrame(columns=col_names)
    
#     var_cluster = all_cluster_acc[f'Cluster {cluster_num}']
    
#     for idx,k in enumerate(var_cluster.keys()):
#         #split to get model, lead, and season
#         split_ = k.split('_')
        
#         df_OUT=df_OUT.append({'Lead':int(split_[0][-1]), \
#                              'Season': split_[-1],'ACC':var_cluster[k]},ignore_index=True)

#     # subset_1 = df_OUT[(df_OUT['Lead'] ==1) | (df_OUT['Lead'] ==2 ) | (df_OUT['Lead'] ==3 )]

#     return(df_OUT)

#%%

all_vals_setup_ETo = {}
for clus_num in np.arange(1,7):
    all_vals_setup_ETo[f'Cluster {clus_num}'] = setup_plot_all_leads(all_cluster_acc=all_cluster_acc_ETo,cluster_num=clus_num)


# r1_r = setup_plot(all_vals['Cluster 1'])


#%%
def plot_lead_week_season_model(all_vals_setup,obs_files,subx_files):
    
    #for each cluster, get all the data prepared for heatmaps
    out_clusters = {}  
    for i in all_vals_setup.keys():
        # all_vals_setup[i]

        cluster_num = int(i.split()[-1])
          
        out_name = return_cluster_name(cluster_num)                
    
        def return_data_for_plot(subset_by_cluster):
            #convert to np array for pcolor mesh
            #season order: Spring, Summer, Fall, Winter
            
            length_of_leads = (all_vals_setup[list(all_cluster_acc_ETo.keys())[0]]['Season'] == 'Winter').sum()
            var_OUT = np.zeros(shape = (4,length_of_leads))
            #Place names in this order for better visual
            var_OUT[0,:] = subset_by_cluster[subset_by_cluster['Season']=='Spring']['ACC']
            var_OUT[1,:] = subset_by_cluster[subset_by_cluster['Season']=='Summer']['ACC']
            var_OUT[2,:] = subset_by_cluster[subset_by_cluster['Season']=='Fall']['ACC']
            var_OUT[3,:] = subset_by_cluster[subset_by_cluster['Season']=='Winter']['ACC']
            
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
                Cols = ['1', '2', '3', '4','3.4']
            else:
                Cols = ['1', '2', '3', '4', '5', '6','3.4','3.5','3.6','4.6']
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
        
#%%
#West region (cluster 1)  
'''min all and max all are for individual days in the first few leads. 
Skill decreases dramatically'''
# RZSM_acc_values_to_plot,max_all_rzsm, min_all_rzsm = plot_lead_week_season_model(all_vals_setup=all_vals_setup_RZSM)
ETo_acc_values_to_plot= plot_lead_week_season_model(all_vals_setup=all_vals_setup_ETo,obs_files=obs_files,subx_files=subx_files)
    
def get_max_and_min_all_clusters(ETo_acc_values_to_plot):
    max_ = 0
    min_ = 1
    
    for i in ETo_acc_values_to_plot.keys():
        if np.max(ETo_acc_values_to_plot[i].values) > max_:
            max_ = np.max(ETo_acc_values_to_plot[i].values)
        if np.min(ETo_acc_values_to_plot[i].values) < min_:
            min_ = np.min(ETo_acc_values_to_plot[i].values)
    return(max_, min_)

max_all, min_all = get_max_and_min_all_clusters(ETo_acc_values_to_plot)     
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
                s.set_title(f'RSMAS CCSM4 \n {var} Anomaly Correlation Coefficient',fontsize=25)
            elif model_NAM1 == 'GMAO':
                s.set_title(f'GMAO GEOS-5 \n {var} Anomaly Correlation Coefficient',fontsize=25)
            elif model_NAM1 == 'ESRL':
                s.set_title(f'ESRL FIMr1p1 \n {var} Anomaly Correlation Coefficient',fontsize=25)
            elif model_NAM1 == 'EMC':
                s.set_title(f'EMC GEFSv12 \n {var} Anomaly Correlation Coefficient',fontsize=25)

    s.set_xlabel('Week Lead',fontsize=25)
    
    plt.savefig(f'{output_season_dir}/all_season_MEM_skill_{var}.tif',dpi=300)
    return(0)

make_save_plots_all_models_seasons_leads(acc_values_to_plot=ETo_acc_values_to_plot,var=var,min_all=min_all,max_all=max_all)








#%% Old code that hasn't been changed to be for the new MEM
# #%%Interannual skill for all grid cells and leads and models 
# # var = 'RZSM'
# def all_skill_by_lead (var, obs_eto_path, obs_rzsm_path,obs_eto_mean,obs_rzsm_mean):
    
#     try:
#         var_yearly = xr.open_dataset(f'{output_nc_dir}/interannual_pearson_skill.nc')
#         return(var_yearly)

#     except FileNotFoundError:
            
#         if var == 'ETo':
#             obs_name = 'ETo_SubX_*.nc4'
#             sub_name = 'ETo_anom'
#             obs_path = obs_eto_path
#         elif var == 'RZSM':
#             obs_name = 'SM_SubX_*.nc4'
#             sub_name = 'RZSM_anom'
#             obs_path = obs_rzsm_path
            
#         #Don't have to load into memory because I'm immediately converting them into a np.array
#         subx_files = xr.open_mfdataset(f'{var}_anomaly_LEAD*.nc4', concat_dim=['S'], combine='nested')
#         obs_files = xr.open_mfdataset(f'{obs_path}/{obs_name}', concat_dim = ['S'], combine = 'nested')

#         # #Get the mean value of all files
#         # subx_mean = xr.open_mfdataset(f'{subdir_for_mean}/{var}_mean*.nc4', concat_dim=['S'], combine='nested') 

#         #Make an empty file to store the final outcomes
#         var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().to_array().to_numpy().squeeze()
#         # var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().copy() #original, but can't s
#         #Because we don't need to slice anything, we can convert it to a numpy array
#         #and use Numba for faster processing
        
#         obs_converted = obs_files.to_array().to_numpy().squeeze() #drop unneeded dimension
#         subx_converted = subx_files.to_array().to_numpy().squeeze() #drop unneeded dimension
        
        
        
#         #Very slow loop (I switched to numba instead)
#          # #Now find pearson correlation by model, lead, and lat/lon
#          # for model in range(var_OUT[f'{sub_name}'].shape[0]):
#          #     print(f'Working on model {model+1} for pearson r correlation')
#          #     for Y in range(var_OUT.Y.shape[0]):
#          #         print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
#          #         for X in range(var_OUT.X.shape[0]):
#          #             for lead in np.arange(7,45,7):
#          #                 var_OUT[f'{sub_name}'][model, lead,Y,X] = \
#          #                     np.corrcoef(subx_files[f'{sub_name}'][:,model, lead, Y, X],obs_files[f'{sub_name}'][:,model, lead, Y, X])[0,1]
                
      
#         def anomaly_correlation_coefficient(var_OUT, subx_converted, obs_converted):
#             #Source ACC:
#             #https://metclim.ucd.ie/wp-content/uploads/2017/07/DeterministicSkillScore.pdf

#             def ACC(FC_anom,OBS_anom):
#                 top = np.nanmean(FC_anom*OBS_anom) #all forecast anomalies * all observation anomalies
#                 bottom = np.sqrt(np.nanmean(FC_anom**2)*np.nanmean(OBS_anom**2)) #variance of forecast anomalies * variance of observation anomalies
#                 ACC = top/bottom
#                 return ACC
            
#         #Now find pearson correlation by model, lead, and lat/lon
#             for model in range(var_OUT.shape[0]):
#                 print(f'Working on model {model+1} for pearson r correlation')
#                 for Y in range(var_OUT.shape[2]):
#                     # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
#                     for X in range(var_OUT.shape[3]):
#                         #anomalies are only present in weekly leads of 7
#                         for lead in np.arange(7,45,7):
                            
#                             '''There is a Zero division error that occurs, to fix this (because numba doesn't like it)
#                             just check and see if the two files have all 0s or np.nans'''
#                             #Find the pearson_r correlation
#                             if len(np.unique(subx_converted[:,model, lead, Y, X])) <=2 or len(np.unique(obs_converted[:,model, lead, Y, X])) <=2:
#                                 break
#                             else:
#                                 var_OUT[model, lead,Y,X] = \
#                                     ACC(subx_converted[:,model, lead, Y, X],obs_converted[:,model, lead, Y, X])
#                             # except ZeroDivisionError:
#                             #     var_OUT[model,lead,Y,X] = np.nan
#             return(var_OUT)
         
        
#         yearly_skill = anomaly_correlation_coefficient(var_OUT, subx_converted, obs_converted)
        
#         #Convert to an xarray object
#         var_yearly = xr.Dataset(
#             data_vars = dict(
#                 Pearson_r_coefficient = (['model','lead','Y','X'], yearly_skill[:,:,:,:]),
#             ),
#             coords = dict(
#                 X = subx_files.X.values,
#                 Y = subx_files.Y.values,
#                 lead = subx_files.lead.values,
#                 model = subx_files.model.values,
        
#             ),
#             attrs = dict(
#                 Description = 'Pearson r correlation over all years and lead times'),
#         ) 
        
#         var_yearly.to_netcdf(f'{output_nc_dir}/interannual_pearson_skill.nc')
        
#     return(var_yearly)



# def plot_interannual_pearsonR(lead_week, mod, var):
#     proj = ccrs.PlateCarree()
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.background_patch.set_facecolor('black') #changes np.nan value colors
#     ax.outline_patch.set_edgecolor('black') #changes the border of the whole plot to black
#     # ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k') #colors ocean blue
    
    
#     ETo = var_yearly.Pearson_r_coefficient.isel(lead=lead_week).sel(model=mod)
#     ETo.plot(
#         transform=proj, subplot_kws={'projection':proj}, cmap='YlOrRd'
#         )
    
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                       linewidth=2, color='gray', alpha=0.5, linestyle='--')
#     gl.xlabels_top = False
#     gl.ylabels_left = False
#     gl.xlines = True #this will remove the longitude lines
#     # gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
    
#     plt.savefig(f'{output_image_dir}/{var}_lead_{lead_week//7}_model{mod}.tif', dpi=300)
#     plt.cla()
#     plt.clf()

# #save plots for each model and each lead week
# for mod in [1,2,3,4]:
#     for lead_week in np.arange(7,49,7):
#         plot_interannual_pearsonR(lead_week,mod,var)