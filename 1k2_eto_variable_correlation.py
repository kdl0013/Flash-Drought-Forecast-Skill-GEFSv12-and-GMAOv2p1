#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TODO: Find the Coefficient of Determination for individual variables 
used to construct reference ET for:
    1.) The entire CONUS
        a.) for each model seperately
        b.) for each weekly lead time
        c.) plot

    2.) For each region individually
        a.) model, weekly lead time, season
@author: kdl
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
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

print(f"PYTHON: {sys.version}")  # PYTHON: 3.8.1 | packaged by conda-forge | (default, Jan 29 2020, 15:06:10) [Clang 9.0.1 ]
print(f" xarray {xr.__version__}")  # xarray 0.14.1
print(f" numpy {np.__version__}")  # numpy 1.17.3
print(f" matplotlib {mpl.__version__}")  # matplotlib 3.1.2


dir1 = 'main_dir'
model_NAM1 = 'model_name'
#%%
# TODO change later
dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
model_NAM1 = 'GMAO'

subX_dir = f'{dir1}/Data/SubX/{model_NAM1}' #where subX model data lies 
os.chdir(f'{subX_dir}')
#We are going to stay in the subx gmao file path


mask_path = f'{dir1}/Data/CONUS_mask/NCA-LDAS_masks_SubX.nc4' #conus mask

#observations
obs_eto_path = f'{dir1}/Data/gridMET/ETo_SubX_values/'


output_image_dir = f'{dir1}/Outputs/ETo_single_variable_skill/{model_NAM1}'
output_season_dir = f'{dir1}/Outputs/coefficient_of_determination/{model_NAM1}/seasonal_skill/'

os.system(f'mkdir -p {output_season_dir}')


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

TODO: Use only the West mask for the West region. Use HP_conus_mask for all other regions.
'''
#%%
#variables for function (test)

# cluster_num = 1
# var='RH'
# var='tasmin'
def all_season_mod_skill(var,cluster_num, mask_path,obs_eto_path):
    '''Split by season, model, lead, then calculate the anomaly correlation 
    coefficient'''
    conus_mask = xr.open_dataset(f'{mask_path}')  
    HP_conus_mask = conus_mask['USDM-HP_mask']
    West_conus_mask = conus_mask['USDM-West_mask']

    obs_name = f'{var}_SubX_*.nc4'
    obs_path = obs_eto_path
        
    if var == 'RH':
        var_subx='RelativeHumidity'
        var_obs = 'dswrf'  #for some reason when I saved the file, it's name is dswrf
    else:
        var_subx=var
            
    #Don't have to load into memory because I'm immediately converting them into a np.array
    subx_files = xr.open_mfdataset(f'{var_subx}_*.nc4', concat_dim=['S'], combine='nested')
    obs_files = xr.open_mfdataset(f'{obs_path}/{obs_name}', concat_dim = ['S'], combine = 'nested')
    
    # subx_files[subx_files == inf]  = np.nan #doesn't work, but should
    
    #Test outputs
    # subx_files.RH[:,0,0,10,10].values
    
    #Only need first 6 leads for subx_files (plus the very first day)
    subx_files = subx_files[f'{var}'][:,:,::7,:,:].to_dataset()
    #Remove the very first day from subx (which isn't actually a lead, just first projection)
    subx_files = subx_files[f'{var}'][:,:,1:,:,:].to_dataset()

    #Lead already removed
    # obs_files = obs_files[f'{var}'][:,:,1:,:,:].to_dataset()
    
    # obs_files = subx_files --only for testing purposes. #TODO: Don't include
    #Change observations to have no infinity (particularly for ETo subx)
    # if var=='RH':
    #     subx_files=subx_files.where(subx_files.RH < -100, np.nan)
    #     subx_files=subx_files.where(subx_files.RH > 1000, np.nan)

    ''''Reassign coordinates because of a random issue where if you don't fix it
    there is one day missing when running the next block of code between subx and obs for each season'''
    obs_files=obs_files.assign_coords(S=subx_files.S.values)
    obs_files = obs_files[f'{var}'][:,:,1:,:,:].to_dataset()
    
    print(f'Units of SubX {var} on same day for single file (it repeats):')
    print(subx_files[f'{var}'][0,3,0,10,10].values)
    print(f'Units of Observations {var} on same day for single file (it repeats):')
    print(obs_files[f'{var}'][0,3,0,10,10].values)
    
    '''PICK SEASON and REGION'''
    #only choose West conus mask file for west cluster
    if cluster_num == 1:
        summer_subx = subx_files[f'{var}'].sel(S=(subx_files['S.season']=='JJA')).where(West_conus_mask == cluster_num)
        fall_subx = subx_files[f'{var}'].sel(S=(subx_files['S.season']=='SON')).where(West_conus_mask == cluster_num)
        winter_subx= subx_files[f'{var}'].sel(S=(subx_files['S.season']=='DJF')).where(West_conus_mask == cluster_num)
        spring_subx= subx_files[f'{var}'].sel(S=(subx_files['S.season']=='MAM')).where(West_conus_mask == cluster_num)
        
        summer_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='JJA')).where(West_conus_mask == cluster_num)
        fall_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='SON')).where(West_conus_mask == cluster_num)
        winter_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='DJF')).where(West_conus_mask == cluster_num)
        spring_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='MAM')).where(West_conus_mask == cluster_num)
        
    elif cluster_num >=2 and cluster_num <=6:
        summer_subx = subx_files[f'{var}'].sel(S=(subx_files['S.season']=='JJA')).where(HP_conus_mask == cluster_num)
        fall_subx = subx_files[f'{var}'].sel(S=(subx_files['S.season']=='SON')).where(HP_conus_mask == cluster_num)
        winter_subx= subx_files[f'{var}'].sel(S=(subx_files['S.season']=='DJF')).where(HP_conus_mask == cluster_num)
        spring_subx= subx_files[f'{var}'].sel(S=(subx_files['S.season']=='MAM')).where(HP_conus_mask == cluster_num)
        
        summer_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='JJA')).where(HP_conus_mask == cluster_num)
        fall_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='SON')).where(HP_conus_mask == cluster_num)
        winter_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='DJF')).where(HP_conus_mask == cluster_num)
        spring_obs = obs_files[f'{var}'].sel(S=(obs_files['S.season']=='MAM')).where(HP_conus_mask == cluster_num)
    
    #Test outputs of RH
    # winter_subx[:,0,::6,10,10].values
    # winter_obs[:,0,0,10,10].values
    # summer_subx[:,0,0,10,10].values

    
    #make lists to iterate through them
    skill_subx = [summer_subx, fall_subx, winter_subx, spring_subx]
    skill_obs = [summer_obs, fall_obs, winter_obs, spring_obs]
    season_name = ['Summer', 'Fall', 'Winter', 'Spring']
    
    #save outputs for each model, lead, season
    output_dictionary = {}
    for season in range(len(skill_subx)):
        if var=='RH' and season == 2:
            '''testing to see why I have np.nan for skill, when there are values in all the datasets
            The reason is becauase there are some infinity values'''
            subx_converted = (skill_subx[season].to_numpy().squeeze())
            subx_converted[subx_converted==inf]=np.nan
        else:  
        #Because we don't need to slice anything, we can convert it to a numpy array
        #and use Numba for faster processing

            subx_converted = (skill_subx[season].to_numpy().squeeze()) #drop unneed dimension
        
        
        obs_converted = (skill_obs[season].to_numpy().squeeze())  #drop unneed dimension
        obs_converted.shape
        #convert to same shape as subx. This meas
        #Make an empty file to store the final outcomes
        var_OUT = np.empty_like(subx_converted)
        var_OUT[:,:,:,:,:] = np.nan
        #Only want to save MME ( lead x Y x X)
        var_OUT = var_OUT[0,0,:,:,:]

        def season_coefficient_of_determination(var_OUT, subx_converted, obs_converted):

            for Y in range(var_OUT.shape[1]):
                # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
                for X in range(var_OUT.shape[2]):
                    #anomalies are only present in weekly leads of 7
                    for lead in range(var_OUT.shape[0]):
                        
                        '''There is a Zero division error that occurs, to fix this (because numba doesn't like it)
                        just check and see if the two files have all 0s or np.nans'''
                        
                        #pearson R from function
                        obs = obs_converted[:,0, lead, Y, X] #Use only a single model because they are all the same
                        subx =np.nanmean(subx_converted[:,:, lead, Y, X],axis=1)
                        # print(np.nanmean(subx))
                        
                        try:
                            var_OUT[lead,Y,X] = pearsonr(obs,subx)[0]**2
                            # r2_score(obs,subx)
                        except ValueError:
                            #Can't perform r2 score on datasets full of np.nan
                            var_OUT[lead,Y,X] = np.nan

            return(var_OUT)
        
        seasonal_skill = season_coefficient_of_determination(var_OUT, subx_converted, obs_converted)
        seasonal_skill=var_OUT
        
        #Now add back to a dictionary for each model/lead week/weason
        
        if var == 'RH':
            lead_range = skill_subx[season].lead.shape[0]
        else:
            lead_range = skill_subx[season].L.shape[0]
        
        for lead_week in range(lead_range):
            seasonal_mod_skill = np.nanmean(seasonal_skill[lead_week,:,:])
            output_dictionary[f'Lead{lead_week+1:02}_{season_name[season]}']= seasonal_mod_skill
            
    return(output_dictionary)
#%%

all_cluster_r2_tasmin = {}
for clus_num in np.arange(1,7):
    print(f'Working on cluster {clus_num} out of 6.')
    all_cluster_r2_tasmin[f'Cluster {clus_num}'] = all_season_mod_skill(var='tasmin',cluster_num=clus_num, mask_path=mask_path,obs_eto_path=obs_eto_path)

#%%
all_cluster_r2_tasmax = {}
for clus_num in np.arange(1,7):
    print(f'Working on cluster {clus_num} out of 6.')
    all_cluster_r2_tasmax[f'Cluster {clus_num}'] = all_season_mod_skill(var='tasmax',cluster_num=clus_num, mask_path=mask_path,obs_eto_path=obs_eto_path)

#%%
all_cluster_r2_RH = {}
for clus_num in np.arange(1,7):
    print(f'Working on cluster {clus_num} out of 6.')
    all_cluster_r2_RH[f'Cluster {clus_num}'] = all_season_mod_skill(var='RH',cluster_num=clus_num, mask_path=mask_path,obs_eto_path=obs_eto_path)

#%%
all_cluster_r2_shortwave_rad = {}
for clus_num in np.arange(1,7):
    print(f'Working on cluster {clus_num} out of 6.')
    all_cluster_r2_shortwave_rad[f'Cluster {clus_num}'] = all_season_mod_skill(var='dswrf',cluster_num=clus_num, mask_path=mask_path,obs_eto_path=obs_eto_path)

#%% Plot anomaly values in pcolormesh

#test
# var_cluster=output_dictionary
# r1_rzsm=output_dictionary

def setup_plot_all_leads(all_cluster_acc,cluster_num):
    #The goal of this function is to take our outputs from anomaly correlation 
    #coefficient, convert them to a dataframe by model, lead, season,
    #to set up for plot
    col_names = ['Lead','Season','Pearson_r2']
    
    df_OUT = pd.DataFrame(columns=col_names)
    
    var_cluster = all_cluster_acc[f'Cluster {cluster_num}']
    
    for idx,k in enumerate(var_cluster.keys()):
        # if idx==1:
        #     break
        #split to get model, lead, and season
        split_ = k.split('_')
        
        df_OUT=df_OUT.append({'Lead':float(split_[0][-2:]), \
                             'Season': split_[-1],'Pearson_r2':var_cluster[k]},ignore_index=True)

    return(df_OUT)


#%%
tasmin_setup = {}
for clus_num in np.arange(1,7):
    tasmin_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(all_cluster_acc=all_cluster_r2_tasmin,cluster_num=clus_num)

tasmax_setup = {}
for clus_num in np.arange(1,7):
    tasmax_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(all_cluster_acc=all_cluster_r2_tasmax,cluster_num=clus_num)

RH_setup = {}
for clus_num in np.arange(1,7):
    RH_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(all_cluster_acc=all_cluster_r2_RH,cluster_num=clus_num)
    
srad_setup = {}
for clus_num in np.arange(1,7):
    srad_setup[f'Cluster {clus_num}'] = setup_plot_all_leads(all_cluster_acc=all_cluster_r2_shortwave_rad,cluster_num=clus_num)
    
#%%
def plot_lead_week_season_model(all_vals_setup,score_name='Pearson_r2'):
    
    def get_max_and_min_all_clusters(all_vals_setup):
        max_ = 0
        min_ = 1
        
        for i in all_vals_setup.keys():
            if max(all_vals_setup[i][f'{score_name}']) > max_:
                max_ = max(all_vals_setup[i][f'{score_name}'])
            if (min(all_vals_setup[i][f'{score_name}'])) < min_:
                min_ = min(all_vals_setup[i][f'{score_name}'])
        return(max_, min_)
    
    max_all, min_all = get_max_and_min_all_clusters(all_vals_setup)     
    
    #for each cluster, get all the data prepared for heatmaps
    out_clusters = {}  
    for i in all_vals_setup.keys():
        # all_vals_setup[i]

        cluster_num = int(i.split()[-1])
          
        out_name = return_cluster_name(cluster_num)                
    
        def return_data_for_plot(subset_by_cluster):
            #convert to np array for pcolor mesh
            #season order: Spring, Summer, Fall, Winter
            var_mod = subset_by_cluster
            
            var_OUT = np.zeros(shape = (4,6))
            #Place names in this order for better visual
            var_OUT[0,:] = var_mod[var_mod['Season']=='Spring'][f'{score_name}']
            var_OUT[1,:] = var_mod[var_mod['Season']=='Summer'][f'{score_name}']
            var_OUT[2,:] = var_mod[var_mod['Season']=='Fall'][f'{score_name}']
            var_OUT[3,:] = var_mod[var_mod['Season']=='Winter'][f'{score_name}']
            
            #Set up for p_color mesh
            x=np.arange(1.5,var_OUT.shape[1]+2)
            x = x-1
            y=np.arange(1.5,var_OUT.shape[0]+2)
            Z=var_OUT[:,:]
            
            Index = ['Spring', 'Summer', 'Fall', 'Winter']
            Cols = ['1', '2', '3', '4', '5', '6',]
            var_OUT = pd.DataFrame(var_OUT,index=Index, columns=Cols)
            
            return(var_OUT,x,y,Z)
        
        #p0-p3 contains the necessary data in a dataframe
        #p0x...p0y etc just contain data that was used for old plotting, not needed anymore
        plot_setup,p0x,p0y,p0z = return_data_for_plot(subset_by_cluster=all_vals_setup[i])


        #add multi_mean to dictionary
        out_clusters[f'{out_name}_MME'] = plot_setup
        
    return(out_clusters,max_all, min_all)
        
#%%
#West region (cluster 1)  
tasmin_acc_values_to_plot,max_all_tasmin, min_all_tasmin = plot_lead_week_season_model(all_vals_setup=tasmin_setup)
tasmax_acc_values_to_plot,max_all_tasmax, min_all_tasmax = plot_lead_week_season_model(all_vals_setup=tasmax_setup)
RH_acc_values_to_plot,max_all_RH, min_all_RH = plot_lead_week_season_model(all_vals_setup=RH_setup)
srad_acc_values_to_plot,max_all_srad, min_all_srad = plot_lead_week_season_model(all_vals_setup=srad_setup)

# max_all_RH=0.75


#%%
def make_save_plots_all_models_seasons_leads(acc_values_to_plot,var,min_all,max_all):
    sorted(acc_values_to_plot)
    
    parameters = {'axes.labelsize': 16,
              'axes.titlesize': 20,
              'xtick.labelsize' : 15,
              'figure.titlesize': 24}
    plt.rcParams.update(parameters)

    # create figure

    cmap = mpl.colors.ListedColormap(plt.cm.YlOrRd(np.linspace(min_all, max_all, 7)))
    # cmap.set_under((.8, .8, .8, 1.0))

    fig, ax = plt.subplots(6, 1, figsize=(20, 10), dpi=300)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    

    region_index= 0
    for idx,region_model in enumerate(sorted(acc_values_to_plot)):
        # if idx == 5:
        #     break
        # print(idx)

        plot_data = acc_values_to_plot[region_model]
        if idx != len(acc_values_to_plot)-1: #only for getting the tick labels to not show up

            #Change region names on the left side
            sns.heatmap(ax=ax[region_index],data = plot_data,
                        cmap=cmap, vmin=min_all, vmax=max_all, annot=True, 
                        fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False,
                        cbar_ax= cbar_ax,xticklabels=False)
            
            ax[idx].set_ylabel(region_model.split('_')[0])
            
            region_index+=1
        else:
            #Change region names on the left side
            sns.heatmap(ax=ax[region_index],data = plot_data,
                        cmap=cmap, vmin=min_all, vmax=max_all, annot=True, 
                        fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False,
                        cbar_ax= cbar_ax,xticklabels=True)
            
            ax[idx].set_ylabel(region_model.split('_')[0])
            ax[idx].set(xlabel='Lead Week')
            
            fig.suptitle(f'gridMET {var} vs. GEOS-5 {var} \nPearson r correlation \n (multi-ensemble mean applied first, then compared with observations)')
            
    return(plt.savefig(f'{output_season_dir}/all_mod_season_lead_{var}_for_ETo_analysis.tif'))


#%%
make_save_plots_all_models_seasons_leads(acc_values_to_plot=tasmin_acc_values_to_plot,var='tasmin',min_all=min_all_tasmin,max_all=max_all_tasmin) 
make_save_plots_all_models_seasons_leads(acc_values_to_plot=tasmax_acc_values_to_plot,var='tasmax',min_all=min_all_tasmax,max_all=max_all_tasmax) 
make_save_plots_all_models_seasons_leads(acc_values_to_plot=RH_acc_values_to_plot,var='RH',min_all=min_all_RH,max_all=max_all_RH) 
make_save_plots_all_models_seasons_leads(acc_values_to_plot=srad_acc_values_to_plot,var='srad',min_all=min_all_srad,max_all=max_all_srad) 





    # fig.tight_layout(rect=[0, 0, .9, 1])
# #%%   
#         count_clus+1
        
#         sns.heatmap(p1,cmap=cmap,cbar=False,ax=ax1, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#         sns.heatmap(p2,cmap=cmap,cbar=False,ax=ax2, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#         sns.heatmap(p3,cmap=cmap,cbar=False,ax=ax3, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
    
    
#     ax0.tick_params(left=False, bottom=False)
#     ax1.tick_params(left=False, bottom=False)
#     ax1.set x
#     ax2.tick_params(left=False, bottom=False)
#     ax3.tick_params(left=False, bottom=False)

#     #%%
#     #Can remove numbers inside of frames
#     fig, (ax0,ax1,ax2,ax3,) = plt.subplots(1, 4, figsize=(29, 10), dpi=300)
#     sns.heatmap(p0,cmap=cmap,cbar=False,ax=ax0, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#     sns.heatmap(p1,cmap=cmap,cbar=False,ax=ax1, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#     sns.heatmap(p2,cmap=cmap,cbar=False,ax=ax2, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#     sns.heatmap(p3,cmap=cmap,cbar=True,ax=ax3, vmin=min_, vmax=max_, annot=True, fmt='.2f', linewidths=2.0, linecolor='black', clip_on=False)
#     ax0.tick_params(left=False, bottom=False)
#     ax1.tick_params(left=False, bottom=False)
#     ax2.tick_params(left=False, bottom=False)
#     ax3.tick_params(left=False, bottom=False)
#     plt.suptitle('test')
#     #%%
#     plt.savefig
    
#     cbar = ax.collections[0].colorbar
#     cbar.set_ticks(np.linspace(1, 7, 2 * ncols + 1)[1::2])
#     cbar.set_ticklabels(range(1, 7 + 1))
        
    

# N = 100
# M = 200
# p = 0.8
# df = pd.DataFrame(np.random.choice([0, 1], (M, N), p=(p, 1 - p)),
#                   columns=sorted((list(range(10)) * N)[0:N]),
#                   index=sorted((list(range(10)) * N)[0:M]))

# ncols = df.index.max() + 1
# cmap = mpl.colors.ListedColormap(plt.cm.turbo(np.linspace(0, 1, ncols)))
# cmap.set_under((.8, .8, .8, 1.0))

# ax = sns.heatmap(df.apply(lambda s: (s.name == s.index) * s * (s.index + 1)), mask=df.eq(0), cmap=cmap, vmin=1)

# cbar = ax.collections[0].colorbar
# cbar.set_ticks(np.linspace(1, ncols, 2 * ncols + 1)[1::2])
# cbar.set_ticklabels(range(1, ncols + 1))
        
#         # create figure
#     fig = plt.figure(figsize=(20, 5), dpi=80)
    
#     # add subplots
#     for i, col in enumerate(iris.columns[:-1], 1):
#         plt.subplot(1, 4, i)
#         ax = sns.boxplot(x='species', y=col, data=iris, hue='species')
#         ax.get_legend().remove()
#         plt.title(col)

# # add legend
# handles, labels = ax.get_legend_handles_labels()
# fig.legend(handles, labels, loc='upper right', ncol=3, bbox_to_anchor=(.75, 0.98))

# # add subtitle
# fig.suptitle('Distribution of floral traits in the species of iris')

    
    
    
    
    
#     fig, (ax0,ax1,ax2,ax3) = plt.subplots(1, 4, figsize=(29, 10), dpi=30)
#     # fig.title('this', loc='center')
#     # plt.title('This is My Legend Title')
#     p0=ax0.pcolormesh(p0x,p0y,p0z,cmap='YlOrRd')
#     ax0.set_title('Model 1')
#     # p0.set_ylabel('this')
    
#     p1=ax1.pcolormesh(p1x,p1y,p1z,cmap='YlOrRd')
#     ax1.set_title('Model 2')
    
#     p2=ax2.pcolormesh(p2x,p2y,p2z,cmap='YlOrRd')
#     ax2.set_title('Model 3')
    
#     p3=ax3.pcolormesh(p3x,p3y,p3z,cmap='YlOrRd')
#     ax3.set_title('Model 4')
    
    
#     fig.colorbar(p0)
    
# #%%    

# sns.heatmap(p0,cmap=plt.cm.YlOrRd,cbar=True)
    
            
#     plt.colorbar()
    
#     c = ax0.pcolor(var_OUT)
#     ax0.set_title('default: no edges')
#  p.axes.set_ylabel('La Niña           Neutral         El Niño')


# plt.pcolormesh(var_OUT)



# x = np.arange(-0.5, 10,5)  # len = 11
# x_len = len(x)
# y = np.arange(4.5, 11, 5)  # len = 7
# y_len = len(y)

# Z = np.random.rand(y_len-1, x_len-1)

# fig, ax = plt.subplots()
# ax.pcolormesh(var_OUT, y, Z)


# Z = np.random.rand(6, 10)


# fig, ax = plt.subplots(2, 1)

# c = ax.pcolor(var_OUT)
# ax.set_title('default: no edges')



# r1_r = r1_r.to_xarray()


# #PLOT
# models = {1: 'Model 1', 2: 'Model 2', 3: 'Model 3', 4: 'Model 4'}

# parameters = {'axes.labelsize': 16,
#           'axes.titlesize': 20,
#           'xtick.labelsize' : 15,
#           'figure.titlesize': 24}

# plt.rcParams.update(parameters)


# # fig, (ax1,ax2) = plt.subplots(3,2)
# # fig=plt.figure()
# p_edfd = r1_r.plot(
#    x='Season', y='Lead', col = 'models',
#    cmap='bwr',
#    cbar_kwargs={'label': ''},
#    figsize=(12,5)
#    )

 
# for i,es in enumerate(enso_states):
#     # ax = fig.add_subplots(2,3,i+1)
#     ax = p_edfd.axes.flat[i]
#     ax.set_title(enso_states.get(i-1, ''))
#     if i ==0:
#         ax.text(-10.08, 4.00,text_,size=30)
#     if i ==1:
#         ax.set_xlabel('')
#         ax.set_title(f'{clus_name}\n\n' +enso_states.get(i-1, ''))
        
#     else:
#         ax.set_xlabel('')
#         ax.set_title(enso_states.get(i-1, ''))
#     ax.set_xticks([])
#     ax.set_yticks([])
# p_edfd.axes[0,0].set_ylabel('SON           JJA            MAM')
# text(0.535, -0.035,f'    <{min_index}%       {min_index}-{max_index}%     >{max_index}%', ha='center', va='center', transform=ax.transAxes,size=13)
# text(-0.51, -0.035,f'    <{min_index}%       {min_index}-{max_index}%     >{max_index}%', ha='center', va='center', transform=ax.transAxes,size=13)
# text(-1.58, -0.035,f'    <{min_index}%       {min_index}-{max_index}%     >{max_index}%', ha='center', va='center', transform=ax.transAxes,size=13)

# #Now plot with pcolor mesh

# np.random.seed(19680801)

# x = np.arange(-0.5, 10,10)  # len = 11
# x_len = len(x)
# y = np.arange(4.5, 11, 5)  # len = 7
# y_len = len(y)

# Z = np.random.rand(y_len-1, x_len-1)

# fig, ax = plt.subplots()
# ax.pcolormesh(x, y, Z)


# #new
# np.random.seed(19680801)
# Z = np.random.rand(6, 10)
# x = np.arange(-0.5, 10, 1)  # len = 11
# y = np.arange(4.5, 11, 1)  # len = 7

# fig, ax = plt.subplots()
# ax.pcolormesh(x, y, Z)
# # subx_files.masked_where(West_conus_mask.West != cluster_num)

# # #Check if files are equivalent (to see if mask worked correclty)
# # np.count_nonzero(np.isnan(subx_files==subx_cluster).to_array())


# # var_yearly=subx_cluster.to_dataset()

# var_yearly = subx_cluster.RZSM_anom[0,0,7,:,:]
# var_yearly=subx_files.RZSM_anom[0,0,7,:,:]

# proj = ccrs.PlateCarree()

# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.background_patch.set_facecolor('white') #changes np.nan value colors
# ax.outline_patch.set_edgecolor('black') #changes the border of the whole plot to black
# # ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k') #colors ocean blue


# # ETo = var_yearly.Pearson_r_coefficient.isel(lead=lead_week).sel(model=mod)
# var_yearly.plot(
#     transform=proj, subplot_kws={'projection':proj}, cmap='YlOrRd'
#     )

# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                   linewidth=2, color='gray', alpha=0.5, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_left = False
# gl.xlines = True #this will remove the longitude lines
# # gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER

# plt.savefig(f'{output_image_dir}/week_lead_{lead_week//7+1}_model{mod}.tif', dpi=300)
# plt.cla()
# plt.clf()



# skill_subx = [summer_subx, fall_subx, winter_subx, spring_subx]
# skill_obs = [summer_obs, fall_obs, winter_obs, spring_obs]
# season_name = ['Summer', 'Fall', 'Winter', 'Spring']

# for season in range(len(skill_subx)):
    
#     #Make an empty file to store the final outcomes
#     var_OUT = skill_subx[season][0,:,:,:,:].to_dataset().to_array().to_numpy().squeeze()

# # var_OUT = subx_files[f'{sub_name}'][0,:,:,:,:].to_dataset().copy() #original, but can't s
# #Because we don't need to slice anything, we can convert it to a numpy array
# #and use Numba for faster processing

#     obs_converted = skill_subx[season].to_numpy().squeeze() #drop unneed dimension
#     subx_converted = skill_obs[season].to_numpy().squeeze() #drop unneed dimension

       
#     @njit
#     def season_skill_assessment(var_OUT, subx_converted, obs_converted):
#     #Now find pearson correlation by model, lead, and lat/lon
#         for model in range(var_OUT.shape[0]):
#             print(f'Working on model {model+1} for pearson r correlation')
#             for Y in range(var_OUT.shape[2]):
#                 # print(f"Working on latitude index {Y} out of {var_OUT.Y.shape[0]}")
#                 for X in range(var_OUT.shape[3]):
#                     for lead in np.arange(7,45,7):
#                         var_OUT[model, lead,Y,X] = \
#                             np.corrcoef(subx_converted[:,model, lead, Y, X],obs_converted[:,model, lead, Y, X])[0,1]
#         return(var_OUT)
    
#     seasonal_skill = season_skill_assessment(var_OUT, subx_converted, obs_converted)
    
#     #Convert to an xarray object
#     var_yearly = xr.Dataset(
#         data_vars = dict(
#             Pearson_r_coefficient = (['model','lead','Y','X'], seasonal_skill[:,:,:,:]),
#         ),
#         coords = dict(
#             X = subx_files.X.values,
#             Y = subx_files.Y.values,
#             lead = subx_files.lead.values,
#             model = subx_files.model.values,
    
#         ),
#         attrs = dict(
#             Description = 'Pearson r correlation for each season and lead times'),
#     ) 
    
#     var_yearly.to_netcdf(f'{output_nc_dir}/{season_name[season]}_pearson_skill.nc')
    
#     return(var_yearly)


# #%%


# xr.where(a.ETo_anom == np.nan, a.ETo_anom, a.ETo_anom) = 0

# '''Test with xskillscore'''
# # a = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SMERGE_SM/SM_SubX_values/SM_SubX_1999-01-10.nc')
# a = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/ETo_anomaly_2012-02-19.nc')
# a.close()
# b = xr.open_dataset('/home/kdl/Insync/OneDrive/NRT_CPC_Internship/Data/SubX/GMAO/ETo_anomaly_2012-02-19.nc')
# b.close()

# rzsm = xr.open_mfdataset('vas*.nc', concat_dim=['S'], combine='nested')
# np.corrcoef(a.ETo_anom[0,0,7,10,10].values,a.ETo_anom[0,0,7,10,10].values)
# # proj = ccrs.PlateCarree()

# only_2_dim = rzsm.vas[::2,0,7,10,10].values
# dim_2 = rzsm.vas[::2,0,7,11,11].values
# len(only_2_dim)
# np.corrcoef(only_2_dim,dim_2)
# np.corrcoef(only_2_dim,dim_2)[0,1]



# # create colormap 
# cmap=plt.get_cmap()

# # use white color to mark 'bad' values
# cmap.set_bad(color='w')
# # gl.xlabel_style = {'size': 15, 'color': 'gray'}
# # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
# # ax.coastlines()

# #example np.corrcoeff
# x_simple = np.array([-2, -1, 0, 1, 2])
# y_simple = np.array([4, 1, 3, 2, 0])
# my_rho = np.corrcoef(x_simple, y_simple)




# np.corrcoef(a.ETo_anom.isel(S=0,model=1,lead=7,Y=10,X=10).values,a.ETo_anom.isel(S=0,model=1,lead=7,Y=10,X=10).values)

# ddd=np.array([1,1,1,1,2])
# ccc = np.array([2,2,2,2,3])
# np.corrcoef(ddd,ccc)

# dataset = a
# # dataset = netcdf_dataset(fname)
# sst = dataset.variables['ETo_anom'][0, 0, 7,:,:]
# lats = dataset.variables['Y'][:]
# lons = dataset.variables['X'][:]

# ax = plt.axes(projection=ccrs.PlateCarree())

# plt.contourf(lons, lats, sst, 60,
#              transform=ccrs.PlateCarree())


# dataset.plot(transform=proj, col='Time', col_wrap=3, robust=True, subplot_kws={'projection':proj})




# import cartopy.crs as ccrs
# import xarray as xr




# dd = xs.pearson_r(a, a, dim='S')
# np.corrcoef(a.to_array(),a.to_array())

# a = xr.DataArray(np.random.rand(5, 3, 3),

#                  dims=['time', 'x', 'y'])

# b = xr.DataArray(np.random.rand(5, 3, 3),

#                  dims=['time', 'x', 'y'])

# xs.pearson_r(a, b, dim='time')

# xs.pearson_r(a,b)

# air1d = c.isel(Y=20, X=10)
# np.nanmean(air1d.to_array())
# air1d.ETo_anom.plot()

# plt.plot(skill['rmm1'])
# plt.title('GMAO-GEOS_V2p1 RMM1 Skill')
# plt.xlabel('Lead Time (Days)')
# plt.ylabel('ACC')

# # dask_gufunc_kwargs.setdefault("allow_rechunk", True)


# HP_conus_mask = xr.open_dataset(HP_conus_path) #open mask file
# West_conus_mask = xr.open_dataset(West_conus_path)

# smerge_file = xr.open_dataset(f'{dir1}/Data/SMERGE_SM/Raw_data/smerge_sm_merged_remap.nc4') #for mask

# '''First, fix RZSM and ETo anomaly SubX files so that the S dimension with no data
# is dropped. This will assist with attempting to do correlation using 
# xr.open_mfdataset'''
# def return_date_list(var):
#     date_list = []
#     for file in sorted(glob(f'{subX_dir}/{var}*.nc4')):
#         date_list.append(file[-14:-4])
#     return(date_list)
        
# init_date_list = return_date_list(var = 'RZSM')    


          

# np.count_nonzero(open_f.where(HP_conus_mask.High_Plains > 0).to_array())
# np.count_nonzero(open_f.where(HP_conus_mask.High_Plains == 0).to_array())  
    

# '''

# USDM-West_mask:
# ---------------

# The U.S. Drought Monitor (USDM) Regions
#  -- For the NCA-LDAS domain

# Data description:

#  -- This dataset includes CONUS-wide (minus Alaska and Hawaii for now)
#     U.S. Drought Monitor regions, delimited by state-political boundaries.
#     These regions were defined by the U.S. Drought Monitor.

# Legend Information:

#  -- Each region has been assigned an integer-based index value here.
#     The corresponding region and integer values include:
#     1:  West region
#     2:  Midwest region
#     3:  HighPlains region
#     4:  South region
#     5:  Southeast region
#     6:  Northeast region

'''
#  -- Note: There are two separate USDM masks - one which should be used
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



#%% EDDI correlation









 