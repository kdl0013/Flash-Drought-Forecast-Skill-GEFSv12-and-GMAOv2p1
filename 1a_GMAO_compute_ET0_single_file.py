#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Notes:
    GMAO calculate reference evapotranspiration (ETo) using ASCE Penman Monteith
    method https://pypi.org/project/refet/ python package : refet
    
    gridMET calculate reference evapotranspiration with the same python package
    
    Calculate EDDI from SubX data (EDDI is already downloaded for historical data)
    


@author: kdl
"""

import xarray as xr
import numpy as np
import os
import datetime as dt
import pandas as pd
from glob import glob
import refet
from multiprocessing import Pool
from scipy.stats import rankdata


dir1 = 'main_dir'
num_processors = int('procs')

# dir1 = '/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
home_dir = f'{dir1}/Data/SubX/GMAO'
os.chdir(home_dir)
#Additional datasets
elevation_dir = f'{dir1}/Data/elevation/'

gridMET_dir = f'{dir1}/Data/gridMET'

smerge_dir = f'{dir1}/Data/SMERGE_SM/Raw_data'
###Files
file_list = os.listdir()

#variables = ['dswrf','tasmax', 'tasmin', 'uas', 'vas']

#TODO: Calculate ETref

#For each date, open each file and compute ETref with et
#All files have the same initialized days (part of the pre-processing that is 
#completed)
def return_date_list(var):
    date_list = []
    for file in sorted(glob(f'{var}*.nc')):
        date_list.append(file[-13:-3])
    return(date_list)
        
init_date_list = return_date_list(var = 'mrso')    
                  #max     #min   #dew   #rad   #wind
#open files [need tasmax, tasmin, tdps, dswrf, windSpeed]

# for _date in init_date_list:
#_date = init_date_list[0]
#%%
'''Compute Reference ET for SubX data'''
def multiProcess_Refet_SubX(_date):
    
    try:
        xr.open_dataset(glob(f'ETo_{_date}.nc4')[0])
        print(f'{_date} already completed for ETo. Saved in {home_dir}.')
    except IndexError:
            
        print(f'Working on date {_date}')
        #Open up each file
        tasmax = xr.open_dataset(glob(f'tasmax*{_date}.nc')[0])
        tasmin = xr.open_dataset(glob(f'tasmin*{_date}.nc')[0])
        tdps = xr.open_dataset(glob(f'tdps*{_date}.nc')[0])
        dswrf = xr.open_dataset(glob(f'dswrf*{_date}.nc')[0])
        windU = xr.open_dataset(glob(f'uas*{_date}.nc')[0])
        windV = xr.open_dataset(glob(f'vas*{_date}.nc')[0])
        elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False)
    
    
    
        def convert_temperature(tasmax, tasmin, dewpoint_temp) -> float:
            '''Convert temperature to Celsius (from Kelvin), calculate vapor pressure (ea)
            from dewpoint temp'''
            tasmax = np.subtract(tasmax, 273.15)
            tasmin = np.subtract(tasmin, 273.15)
            dewpoint_temp = np.subtract(dewpoint_temp, 273.15)
            ea = 0.6112 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
            
            return(tasmax, tasmin, dewpoint_temp, ea)
        
        tasmax, tasmin, dewpoint_temp, ea = convert_temperature(tasmax = tasmax, tasmin = tasmin,
                                                         dewpoint_temp = tdps)     
    
        def convert_dswrf_RN(shortwave_radiation) -> float:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            shortwave_rad = shortwave_radiation / 1000000
            return(shortwave_rad)
    
        dswrf = convert_dswrf_RN(shortwave_radiation = dswrf)
    
        def compute_windSpeed(windU, windV) -> float:
            ''' Returns average daily wind at 10-m height m/s'''
            #Square u and v and take the square root to get windspeed
            windSpeed = xr.ufuncs.sqrt((np.square(windU.uas) + np.square(windV.vas)))
            return (windSpeed)
    
        windSpeed = compute_windSpeed(windU = windU, windV = windV)
        
        output_f = xr.zeros_like(tasmax)
        
        for i_mod in range(tasmax.tasmax.shape[1]):
            for i_lead in range(tasmax.tasmax.shape[2]):
                #changed these 2 lines of code from under i_X
  
                     
                for i_Y in range(tasmax.tasmax.shape[3]):
                    for i_X in range(tasmax.tasmax.shape[4]):
                        # print(i_X)
                        date_val = pd.to_datetime(tasmax.S.values[0]) + dt.timedelta(days=i_lead)
                        day_doy = date_val.timetuple().tm_yday #julian day
                    
                        output_f.tasmax[0,i_mod, i_lead, i_Y, i_X] = \
                            (refet.Daily(tmin = tasmin.tasmin[0,i_mod, i_lead, i_Y, i_X].values, \
                                        tmax = tasmax.tasmax[0,i_mod, i_lead, i_Y, i_X].values, \
                                        ea = ea.tdps[0,i_mod, i_lead, i_Y, i_X].values, \
                                        rs = dswrf.dswrf[0,i_mod, i_lead, i_Y, i_X].values, \
                                        uz = windSpeed[0,i_mod, i_lead, i_Y, i_X].values, \
                                        zw = 10, elev = elevation.data[0,i_Y, i_X].values,
                                        lat = tasmax.Y[i_Y].values,
                                        doy = day_doy, method = 'asce',
                                        input_units={'tmin': 'C', 'tmax': 'C', \
                                                     'rs': 'Langleys', 'uz': 'm/s', \
                                                         'lat': 'deg'}).etr())[0]                            
    
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo = (['S', 'model','lead','Y','X'], output_f.tasmax.values),
            ),
            coords = dict(
                X = output_f.X.values,
                Y = output_f.Y.values,
                lead = output_f.L.values,
                model = output_f.M.values,
                S = output_f.S.values
            ),
            attrs = dict(
                Description = 'Reference crop evapotranspiration (mm/day)'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{home_dir}/ETo_{_date}.nc4', mode ='w')
        print(f'Saved file into {home_dir}.')
#%%    

if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_Refet_SubX, init_date_list)


#%%
'''Compute Reference ET for gridMET data.
 Because gridMET refET in its orginial state only contains 2 significant figures,
 i.e., 1.3 mm/d, we need to compute refET using gridMET data. But also, all gridmet
 data is 2 sig. fig. length, but this will create a float value.
 '''
def Refet_gridMET():
    
    try:
        xr.open_dataset(f'{gridMET_dir}/ETo_gridMET_merged.nc')
        print('Already completed gridMET Reference ET.')
    except FileNotFoundError:
        
        print('Working on computing ET ref from gridMET variables.')
        slice1 = "1999-01-10"
        slice2 = "2016-02-09"
        # gridMET_file = gridMET_file.sel(day=slice(f"{slice1}",f"{slice2}"))
        #Open up each file
        tasmax = xr.open_dataset(f'{gridMET_dir}/tmmx_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        tasmin = xr.open_dataset(f'{gridMET_dir}/tmmn_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        tavg = xr.open_dataset(f'{gridMET_dir}/tavg_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #Kelvin
        dswrf = xr.open_dataset(f'{gridMET_dir}/srad_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #already in W/m2
        wind = xr.open_dataset(f'{gridMET_dir}/vs_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
        RHmin = xr.open_dataset(f'{gridMET_dir}/rmax_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
        RHmax = xr.open_dataset(f'{gridMET_dir}/rmin_remap_final.nc').sel(day=slice(f"{slice1}",f"{slice2}")) #m/s2, 10 ft.
        # gridMET_file = gridMET_file.sel(day=slice("1999-01-11","2016-02-09"))
            
        elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False)
        
        
        def convert_temperature(tasmax, tasmin) -> float:
            '''Convert temperature to Celsius (from Kelvin), calculate vapor pressure (ea)
            from dewpoint temp'''
            tasmax = np.subtract(tasmax, 273.15)
            tasmin = np.subtract(tasmin, 273.15)
            
            return(tasmax, tasmin)
        #Call function for temperature files
        tasmax, tasmin = convert_temperature(tasmax = tasmax, tasmin = tasmin)
                
        test_temp = 20
        def esat_vapor_pressure(tavg, test_temp) -> float:
            '''Tetens Formula: Need temp in Celsius. Returns values in kPa
            Test temp should return 2339Pa or 2.34kPa'''
            esat = 0.6112 * np.exp((17.27*tavg)/(237.3 + tavg))
            esat_test = 0.6112 * np.exp((17.27*test_temp)/(237.3 + test_temp))
            delta = (2504 * esat) / (tavg + 237.3)**2
            return(esat, delta, esat_test)
        
        #Calculate saturated vapor pressure:
        esat, delta, esat_test = esat_vapor_pressure(tavg, test_temp)
        #Calculate vapor pressure:
        min_out, DN, DN = esat_vapor_pressure(tasmin, test_temp)
        max_out ,DN, DN = esat_vapor_pressure(tasmax, test_temp)
        
        ea = xr.zeros_like(tavg)
        ea = ((min_out.air_temperature * (RHmax.relative_humidity/100)) + (max_out.air_temperature * (RHmin.relative_humidity/100))) / 2
            

        def convert_dswrf_RN(shortwave_radiation) -> float:
            '''Returns shorwave radiation in MJ/m2 from W/m2'''
            shortwave_rad = shortwave_radiation / 1000000
            return(shortwave_rad)
    
        dswrf = convert_dswrf_RN(shortwave_radiation = dswrf)
    
        output_f = xr.zeros_like(tasmax)

        for day_ in range(tasmax.day.shape[0]):
            for lat in range(tasmax.lat.shape[0]):
                for lon in range(tasmax.lon.shape[0]):
                    day_doy = pd.to_datetime(tasmax.day[day_].values).timetuple().tm_yday #julian day
                    
                    output_f.air_temperature[day_, lat, lon] = \
                        (refet.Daily(tmin = tasmin.air_temperature[day_, lat, lon].values, \
                                    tmax = tasmax.air_temperature[day_, lat, lon].values, \
                                    ea = ea[day_, lat, lon].values, \
                                    rs = dswrf.surface_downwelling_shortwave_flux_in_air[day_, lat, lon].values, \
                                    uz = wind.wind_speed[day_, lat, lon].values, \
                                    zw = 10, elev = elevation.data[0,lat, lon].values,
                                    lat = tasmax.lat[lat].values,
                                    doy = day_doy, method = 'asce',
                                    input_units={'tmin': 'C', 'tmax': 'C', \
                                                 'rs': 'Langleys', 'uz': 'm/s', \
                                                     'lat': 'deg'}).etr())[0]                            
    
        #Convert to an xarray object
        var_OUT = xr.Dataset(
            data_vars = dict(
                ETo_gridmet = (['day', 'lat','lon'], output_f.air_temperature.values),
            ),
            coords = dict(
                day = output_f.day,
                lat = output_f.lat,
                lon = output_f.lon
            ),
            attrs = dict(
                Description = 'Reference crop evapotranspiration (mm/day)'),
        )                    
        
        #Save as a netcdf for later processing
        var_OUT.to_netcdf(path = f'{gridMET_dir}/ETo_gridMET_merged.nc', mode ='w')
        print(f'Saved file into {gridMET_dir}.')
#%%    

#Run function, saved
Refet_gridMET()
#%% Compute EDDI for SubX

'''Steps for EDDI calculation. 
1.) For each date:
    2.) Go thorugh all files and find the 1 week summed ETo for the same dates:
        2a.) Only append the first file with rank and EDDI
        3.) Then rank them according to which weeks have the highest summation:
        4.) Calculate Tukey plotting position
        5.) Calculate EDDI
    
    
To make the multiprocessing function faster (because it appears that accessing the same file is
actually slowing down the processing), I am creating many instances of the 
ETo_ SubX files into new directories. Then I can np.random a number and then
we can pull from that directory.
    
'''

# for _date in init_date_list:   
def multiProcess_EDDI_SubX(_date):
 
    try:
        xr.open_dataset(glob(f'{home_dir}/EDDI_{_date}_completed.nc4')[0])
        print(f'{_date} already completed for EDDI. Saved in {home_dir}.')
    except IndexError:
        
        #Testing different seed for multiprocessing
        np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
        
        #to make a 1-week EDDI, we will need to find the sum of 7 days. This
        #will assist with that <week_lead>
        week_lead = 6
        
        random_directory = np.random.choice(range(1,41))
        print(f'Picked random directory ETo{random_directory}')
        
        os.chdir(f'{home_dir}/ETo{random_directory}')
        
        print(f'Working on date {_date} for EDDI and saving into {home_dir}')
        #Open file, get dates
        Et_ref = xr.open_dataset(f'ETo_{_date}.nc4')
        
        #Make a copy to append to
        Et_out = xr.zeros_like(Et_ref)
        
        #file for SMERGE to not have to run needless loops within CONUS
        smerge_file = xr.open_dataset(f'{smerge_dir}/smerge_sm_merged_remap.nc4')
      
        #Now we have a file with date 1999-01-11. Now we need to find all files
        #with dates that are within 7 days of the current date val
        #get the dates into a list
        a_date_in= Et_ref.ETo.lead.values
        #get the start date
        a_start_date = pd.to_datetime(Et_ref.S.values[0])
        #Add the dates based on index value in date_in
        a_date_out = []
        
        for a_i in range(len(a_date_in)):
            a_date_out.append(a_start_date + dt.timedelta(days = a_i))
            
        start_julian = pd.to_datetime(Et_ref.S[0].values).timetuple().tm_yday #julian day
                                    
        a_julian_out = [start_julian + i for i in range(len(a_date_out))]
        
        INdate_for_month = dt.datetime(int(_date[0:4]),int(_date[5:7]),int(_date[8:10]))

        #Now convert to julian date and append coordinates
        Et_ref2 = Et_ref.assign_coords(lead = a_julian_out)
        
        for model in range(Et_ref.ETo.shape[1]):
            print(f'Working on model {model}')
            for i_Y in range(Et_ref.ETo.shape[3]):
                for i_X in range(Et_ref.ETo.shape[4]):

                    #Don't loop over the values that shouldn't have values.  We will have 4 missing grid cells,
                    #but the computational cost is well worth it. SMERGE indexes 38,6 and 38,7 are EDDI values that we can compute within CONUS
                    '''After calculating EDDI, the researcher discovered that this next line of code added some grid cell calculations
                    that are outside of CONUS. Because this was so computationally taxing, I will just add a new script to change all of the
                    values that are outside of CONUS grid area.'''
                    if ((np.count_nonzero(np.isnan(smerge_file.RZSM[0,i_Y,i_X].values)) ==1) and ((i_X != 38 and i_Y != 6) or (i_X !=38 and i_Y !=7))):
                        Et_out.ETo[0,model,:,i_Y,i_X] = np.nan
                         
                    else:
                                                
                        #append 7-day summation from all files to a new dictionary
                        summation_ETo = {}
                    
                        for val in a_julian_out:
                            #You must val + week_lead because with EDDI you need 7-day looking backwards into the past.
                            #We can

                            if (len(Et_ref2.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=model, X=i_X, Y=i_Y).values)) < 7:
                                break
                            else:
                                summation_ETo[f'{week_lead+val}']=[]
                                summation_ETo[f'{week_lead+val}'].append({f'{_date}':np.nansum(Et_ref2.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=model, X=i_X, Y=i_Y).values)})   
                                # print(len(Et_ref2.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=model, X=i_X, Y=i_Y).values))
                        '''Now we have created a dictionary that contains the 7-day sum of Eto by index:
                            Next we will append to each julian day value in the key in the dictionary with
                            the same julian day from all files'''
                            
                            
                            #TODO: Need to add another dictionary that keeps up with the initialized date of the file
                            #this will assist with not having to re run 1000 files each time.
                            
                        
                        for file in sorted(glob('ETo*.nc4')):
                            #Remove unnecessary files that won't contain any of the dates, if difference in months
                            #is >=2, then skip that file because of 45 day lead time.
                            
                            OUTdate_for_month = dt.datetime(int(file[4:8]),int(file[9:11]),int(file[12:14]))
                                                        
                            num_months = (OUTdate_for_month.month - INdate_for_month.month)
                            
                            if num_months <3: #skip months that aren't within 3 months
                                #Dont' re-open the same file
                                if file[-14:-4] != _date:
                                    
                                    open_f = xr.open_dataset(file)
                                    
                                    '''Convert lead dates to a vector, then add it back into a netcdf
                                    because I cannot convert np.datetime to a pd.datetime, I will need
                                    to iterate over each of the dates in a list from Et_ref'''
                                    
                                    #get the dates into a list
                                    date_in= open_f.ETo.lead.values
                                    #get the start date
                                    start_date = pd.to_datetime(open_f.S.values[0])
                                    #Add the dates based on index value in date_in
                                    date_out = []
                                    for i in range(len(date_in)):
                                        date_out.append(start_date + dt.timedelta(days = i))
                                    
                                    #Convert to julian date
                                    end_julian = pd.to_datetime(open_f.S[0].values).timetuple().tm_yday #julian day
                                    
                                    b_julian_out = [end_julian + i for i in range(len(date_out))]
    
                                    Et_ref_open_f = open_f.assign_coords(lead = b_julian_out)
                                    
                                    '''Now we need to append to the dictionary with the same julian date values'''
                                    for val in b_julian_out:
                                        try:
                                            summation_ETo[f'{week_lead+val}'].append({f'{file[-14:-4]}':np.nansum(Et_ref_open_f.ETo.sel(lead=slice(val,week_lead+val)).isel(S=0, model=model, X=i_X, Y=i_Y).values)})   
                                        except KeyError:
                                            pass
                                        
                                    open_f.close()
                                '''Now the dictionary is filled with all values from all files, we need to sort
                                each julian day by having the largest values as number 1 ranking. Then we can append the 
                                Et_out file with the proper julian day and the EDDI value'''
                             
                        '''Now that we have the total number of files which have all possible days, 
                        we need to rank them'''
                        #chunk start
                        out_eddi_dictionary = {}
                        for idx,julian_date in enumerate(summation_ETo):
                            
                            out_eddi_dictionary[f'{julian_date}']=[]
                            
                            subset_by_date = summation_ETo[f'{julian_date}']
                            
                            #When looking at each julian date, now create a ranked list
                            list_rank = []
                            for date_init in range(len(subset_by_date)):
                                list_rank.append(list(subset_by_date[date_init].values())[0])
                            
                            
                                ranking = rankdata([-1 * i for i in list_rank]).astype(int)
                                tukey = (ranking - 0.33)/ (len(ranking) + 0.33)
                                
                            out_eddi = []
                            '''Now we can calculate EDDI'''
                            for valu in tukey:
                                if valu <= 0.5:
                                    w = out_eddi.append(np.sqrt(-2*np.log(valu)))
                                else:
                                    w = out_eddi.append(1-valu)  
                                        
                            #constants
                            c0 = 2.515517
                            c1 = 0.802853
                            c2 = 0.010328
                            d1 = 1.432788
                            d2 = 0.189269
                            d3 = 0.001308
                                
                            final_out_eddi = []
                            for w in out_eddi:
                                eddi = final_out_eddi.append(w - ((c0 + c1*w + c2*w**2)/(1 + d1*w + d2*w**2 + d3*w**3)))
                             
                            out_eddi_dictionary[f'{julian_date}'].append(final_out_eddi)
                                #Make a new dictionary instead of trying to append to a new one
                       #chunk end     
                        '''Instead of re-looping through all of the files (very slow), we can start appending to 
                        files one by one with the data that we have already collected.
                        
                        The above chunk calculated EDDI, next step is to append to the summation_ETo file becuase that file
                        has the initialized dates for each EDDI julian day summation'''
                        
                        #TODO: Replace summation_ETo values with EDDI values (much much harder than it sounds)
                        
                        final_out_dictionary_all_eddi = {}
                        
                        for idx,julian_dattt in enumerate(summation_ETo):
                            final_out_dictionary_all_eddi[f'{julian_dattt}'] = [] #create an empty list to append to
                            sub_list = summation_ETo[f'{julian_dattt}'] #choose only the summation ETo with correct julian date
                            sub_keys = [] #initialize a list to keep up with the correct julian date and the actual init dates (because each init date varies with number of samples)

                            #sub list contains a dictionary for each julian date in current loop and the values of ETo
                            #Save the init date values for each julian date
                            for idxxx, dictonary in enumerate(sub_list):
                                sub_keys.append({list(dictonary.keys())[0] :out_eddi_dictionary[f'{julian_dattt}'][0][idxxx]})
                            
                            
                            final_out_dictionary_all_eddi[f'{julian_dattt}'].append(sub_keys) 
                        
                        
                        
                        '''Now that we have final_out_dictionary_all_eddi which contains the specific values for each init date for the currenly looped X,Y grid cell, 
                        we can append to a new file (or an already existing file)'''
                        
                        #I think I only need to do 1 year's worth of data at most, make all files first (they are empty)
                        for file_2 in init_date_list:
                            
                            
                            #make files
                            #need to convert the initialized date into julian date and add 45
                            #init date is actually the date-1 (so 01-10-1999 only has data beginning 01-11-1999 with 44 more days)
                            
                            st_day =  pd.to_datetime(file_2)
                            day_doyey = st_day.timetuple().tm_yday + 1 #julian day

                            #add more julian dates
                            day_julian_a = [pd.to_datetime(day_doyey) + dt.timedelta(days=i) for i in range(45)]
                            day_julian_b = [i.timetuple().tm_yday for i in day_julian_a]                            

                            
                            S_values = [pd.to_datetime(file_2)+ dt.timedelta(days=1),pd.to_datetime(file_2)]
                            
                            
                            var_OUT = xr.Dataset(
                            data_vars = dict(
                                    EDDI = (['S', 'model','lead','Y','X'], Et_out.ETo.values),
                                ),
                                    coords = dict(
                                    X = Et_ref.X.values,
                                    Y = Et_ref.Y.values,
                                    lead = day_julian_b,
                                    model = Et_ref.model.values,
                                    S = S_values
                                ),
                                attrs = dict(
                                    Description = 'Evaporative demand drought index from SubX leads'),
                            )                    
                            
                            #Save as a netcdf for later processing
                            var_OUT.to_netcdf(path = f'{home_dir}/EDDI_{file_2}.nc4', mode ='w')
                              

                            '''Now that we have created new files, we can append each file with the data that was found'''
                            for i_val in final_out_dictionary_all_eddi:
                                #for some reason it's a list in a list, this fixes that
                                BB = final_out_dictionary_all_eddi[f'{i_val}'][0]
                                
                                for dic_values in BB:
                                    #Open up the file
                                    eddi_open = xr.open_dataset(f'{home_dir}/EDDI_{list(dic_values.keys())[0]}.nc4')
                                    
                                    print(dic_values)
                                
                                
                                #Open each init date_file and append to the proper date
                                
                                
                                
                                
                                #Find the index of the i_val
                                lead_values = eddi_open.lead.values 
                                index_val=np.where(lead_values == int(i_val))[0][0]

                                
                                eddi_open.EDDI[0,model,index_val,i_Y,i_X] = 



                                                             
                            try:
                                #try to append to an existing file
                                file_open = xr.open_dataset(f'{home_dir}/EDDI_{_date}.nc4', mode = 'w')
                                
                                    # for ii in range(len(final_out_dictionary_all_eddi)):
                                    
                                   
                                    
                                    
                                    
                                    
                                    
                                    #
                                    for b in BB:
                                        print(b)
                                        
                                    
                                    
                                    print(i)
                                for 
                            
                            #TODO: add code to append eddi values
                            
                            
                        except IndexError:
                            #create new file because it doesn't exist to an xarray object
                            var_OUT = xr.Dataset(
                                data_vars = dict(
                                    EDDI = (['S', 'model','lead','Y','X'], Et_out.ETo.values),
                                ),
                                coords = dict(
                                    X = Et_ref.X.values,
                                    Y = Et_ref.Y.values,
                                    lead = a_julian_out,
                                    model = Et_ref.model.values,
                                    S = Et_ref.S.values
                                ),
                                attrs = dict(
                                    Description = 'Evaporative demand drought index from SubX leads'),
                            )                    
                            
                            
                            
                            
                            #Now append any and all possible values into this file, only use the 
                            
                            
                            #Save as a netcdf for later processing
                            var_OUT.to_netcdf(path = f'{home_dir}/EDDI_{_date}.nc4', mode ='w')
                            print(f'Saved _date into {home_dir}.')       
                                                                        
                       
                            
                            for idxxx, dictonary in enumerate(sub_list):
                                final_out_dictionary_all_eddi[f'{julian_dattt}'].append({sub_keys[idxxx]:out_eddi_dictionary[f'{julian_dattt}'][0][idxxx]})   
                                
                                
                            #Now we have the dates, now combine the julian date/init_date from sub_keys/sub_list to  out_eddi_dict
                            
                            
                            
                            
                            
                            
                            print(idx,val)
                        
                        
                                    
                  
                            
                            '''Now append to an empty file'''
                            Et_out.ETo[0,model,week_lead+idx,i_Y,i_X] = eddi
                            
                            
                    
       
             
                        
#%%        
#TODO: see below        
'''For multiprocessing test, need to move the first file created into a new folder.
#Then run multiprocessing to see if it will re-create the same file. The rationale is that                 
#multiprocessing will be pulling from the same file and I'm not sure if multiprocessing is allowed
todo that'''
if __name__ == '__main__':
    p = Pool(num_processors)
    p.map(multiProcess_EDDI_SubX, init_date_list)
    
    #%%
                   

'''This old code was a different formulation for reference ET. It didn't produce
meaningful results across CONUS (despite needing data more data as input than the current). 
This is the ASCE reference ET package'''
#  -- don't run    

# soil_moisture = xr.open_dataset(glob('mrso*')[0]).stack(grid = ['X','Y'])
# #Find nan locations for later use
# nan_dataset = np.where(np.isnan(soil_moisture.READ_DESCRIPTION_ATTRS.values), np.nan, 0)

# #Open file and change all values of 0.0 to nan
# temperature_1 = (xr.open_dataset(glob('tas*')[0]).stack(grid = ['X','Y']))
# temperature_1 = temperature_1.where(temperature_1 != 0)

# windV = (xr.open_dataset(glob('vas*')[0]).stack(grid = ['X','Y']))
# windV = windV.where(windV != 0)

# windU = (xr.open_dataset(glob('uas*')[0]).stack(grid = ['X','Y']))
# windU = windU.where(windU != 0)

# precip = (xr.open_dataset(glob('pr*')[0]).stack(grid = ['X','Y']))
# precip = precip.where(precip != 0)

# dewpoint_temp_1 = (xr.open_dataset(glob('tdps*')[0]).stack(grid = ['X','Y']))
# dewpoint_temp_1= dewpoint_temp_1.where(dewpoint_temp_1 != 0)

# shortwave_radiation = (xr.open_dataset(glob('dswrf*')[0]).stack(grid = ['X','Y']))
# shortwave_radiation = shortwave_radiation.where(shortwave_radiation != 0)

# elevation = xr.open_dataset(f'{elevation_dir}/elev_regrid.nc', decode_times=False).stack(grid = ['lon','lat'])
# elevation = elevation.data[0,:]

# dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values
# #%%Temp is good, esat is good, delta looks good, 
# #Step 1, Compute ET0 using the lagged average ensemble

# #Constants
# Cn = 1600 #mm/d
# Cd = 0.38 #mm/s

# def convert_temperature(temperature, dewpoint_temp) -> float:
#     '''Convert temperature to Celsius (from Kelvin)'''
#     temperature = np.subtract(temperature_1, 273.15)
#     dewpoint_temp = np.subtract(dewpoint_temp_1, 273.15)
    
#     return(temperature, dewpoint_temp)

# temperature, dewpoint_temp = convert_temperature(temperature = temperature_1, \
#                                                  dewpoint_temp = dewpoint_temp_1)

# # #Check outputs
# # print(temperature.READ_DESCRIPTION_ATTRS.values)
# # print(temperature.READ_DESCRIPTION_ATTRS[0,200,0].values)
# # print(np.nanmax(temperature.READ_DESCRIPTION_ATTRS.values))
# # print(np.nanmin(temperature.READ_DESCRIPTION_ATTRS.values))

# # print(dewpoint_temp_1.READ_DESCRIPTION_ATTRS.values)
# # print(dewpoint_temp.READ_DESCRIPTION_ATTRS[0,200,0].values)
# # print(np.nanmax(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
# # print(np.nanmin(dewpoint_temp.READ_DESCRIPTION_ATTRS.values))
# ##############################################################################
# test_temp = 20
# def esat_vapor_pressure(temperature, test_temp) -> float:
#     '''Tetens Formula: Need temp in Celsius. Returns values in kPa
#     Test temp should return 2339Pa or 2.34kPa'''
#     esat = 0.6112 * np.exp((17.27*temperature)/(237.3 + temperature))
#     esat_test = 0.6112 * np.exp((17.27*test_temp)/(237.3 + test_temp))
#     delta = (2504 * esat) / (temperature + 237.3)**2
#     return(esat, delta, esat_test)

# esat, delta, esat_test = esat_vapor_pressure(temperature = temperature, test_temp=test_temp)

# print()
# print(f'Esat_test kPa is {esat_test} at {test_temp} C.')
# print()
# # print(esat.READ_DESCRIPTION_ATTRS.values)
# print(f'Max esat value is {np.nanmax(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min esat value is {np.nanmin(esat.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(delta.READ_DESCRIPTION_ATTRS.values)
# print(f'Max delta value is {np.nanmax(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min delta value is {np.nanmin(delta.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# ##############################################################################
# test_dewpoint_temp = 16
# def e_vapor_pressure(esat, dewpoint_temp, test_dewpoint_temp) -> float:
#     '''Need temp in Celsius. Returns vapor pressure in Pa and humidity
#     values in percent (e.g., 0.78)
#     ea_test should return 1.81kPa'''
#     ea = 0.6108 * np.exp((17.27*dewpoint_temp)/(237.3 + dewpoint_temp))
#     ea_test = 0.6108 * np.exp((17.27*test_dewpoint_temp)/(237.3 + test_dewpoint_temp))
#     rel_humidity = np.divide(ea,esat)
#     return(ea, rel_humidity, ea_test)

# ea, rel_humidity, ea_test = e_vapor_pressure(esat = esat, dewpoint_temp = dewpoint_temp,test_dewpoint_temp = 16)

# print(f'Ea_test kPa is {ea_test} at {test_dewpoint_temp} C.')
# print()
# # print(ea.READ_DESCRIPTION_ATTRS.values)
# print(f'Max ea value is {np.nanmax(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min ea value is {np.nanmin(ea.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(rel_humidity.READ_DESCRIPTION_ATTRS.values)
# print(f'Max rel_humidity value is {np.nanmax(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
# print()
# print(f'Min rel_humidity value is {np.nanmin(rel_humidity.READ_DESCRIPTION_ATTRS.values)}%')
# print()
# ##############################################################################
# def convert_dswrf_RN(shortwave_radiation) -> float:
#     '''Returns shorwave radiation in MJ/m2 from W/m2'''
#     shortwave_rad = shortwave_radiation / 1000000
#     return(shortwave_rad)

# shortwave_rad = convert_dswrf_RN(shortwave_radiation = shortwave_radiation)

# # print(shortwave_rad.READ_DESCRIPTION_ATTRS.values)
# print(f'Max shortwave_rad value is {np.nanmax(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
# print()
# print(f'Min shortwave_rad value is {np.nanmin(shortwave_rad.READ_DESCRIPTION_ATTRS.values)} MJ/m2')
# print()
# ##############################################################################
# #Square both u and v and then take the square root
# def compute_windSpeed(windU, windV) -> float:
#     ''' Returns average daily wind at 10-m height m/s'''
#     #Square u and v and take the square root to get windspeed
#     windSpeed = xr.ufuncs.sqrt(np.square(windU) + np.square(windV))
#     return (windSpeed)

# windSpeed = compute_windSpeed(windU = windU, windV = windV)

# # print(windSpeed.READ_DESCRIPTION_ATTRS.values)
# print(f'Max windSpeed value is {np.nanmax(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
# print()
# print(f'Min windSpeed value is {np.nanmin(windSpeed.READ_DESCRIPTION_ATTRS.values)} m/s')
# print()


# ##############################################################################
# def mean_atmospheric_pressure(elevation) -> float:
#     pressure = xr.zeros_like(temperature)

#     for k in range(elevation.data.shape[0]):
#         if k == 0 and i == 0 and j == 0:
#             print(f'Pressure at elevation {elevation[k].values} meters is {101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)} kPa.')
#         pressure.READ_DESCRIPTION_ATTRS[:,:,k] = 101.3 * (((293 - 0.0065*elevation[k].values)/293)**5.26)
#         # pressure = pressure.where(temperature_1 != 0)

#     psychro = 0.000665 * pressure
#     '''Returns atmospheric pressure in kPa. And returns psychrometric constant'''
    
#     return(pressure, psychro)

# pressure, psychro = mean_atmospheric_pressure(elevation = elevation)


# # print(pressure.READ_DESCRIPTION_ATTRS.values)
# print()
# print(f'Max pressure value is {np.nanmax(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()
# print(f'Min pressure value is {np.nanmin(pressure.READ_DESCRIPTION_ATTRS.values)} kPa')
# print()

# # print(psychro.READ_DESCRIPTION_ATTRS.values)
# print(f'Max psychro value is {np.nanmax(psychro.READ_DESCRIPTION_ATTRS.values)}')
# print()
# print(f'Min psychro value is {np.nanmin(psychro.READ_DESCRIPTION_ATTRS.values)}')
# print()

# #%%
# #TODO
# '''NEED to convert all negative numbers to 0'''
# def compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, \
#                      esat, ea):
    
#     # ETref = xr.zeros_like(windSpeed)
#     '''Returns Reference Evapotranspiration in mm/d'''
#     ETref = ((0.408 * delta * shortwave_radiation) + (psychro*(Cn/temperature+273) * windSpeed * \
#         (esat - ea))) / (delta + (psychro * (1 + Cd * windSpeed)))
        
#     return(ETref)
    
# ETref_func = compute_ETo_crop(delta, shortwave_radiation, psychro, Cn, Cd, temperature, windSpeed, esat, ea)

# ETref_func.READ_DESCRIPTION_ATTRS.values

# # print(ETref_func.READ_DESCRIPTION_ATTRS.values)
# print()
# print(f'Max ETref_func value is {np.nanmax(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
# print()
# print(f'Min ETref_func value is {np.nanmin(ETref_func.READ_DESCRIPTION_ATTRS.values)}')
# print()


# ETref_unstack = ETref_func.unstack('grid').transpose('realization', 'time_series_dates', 'Y', 'X')

# #Convert to an xarray object
# var_OUT = xr.Dataset(
#     data_vars = dict(
#         ETo = (['realization','time_series_dates','Y','X'], ETref_unstack.READ_DESCRIPTION_ATTRS.values),
#     ),
#     coords = dict(
#         X = ETref_unstack.X.values,
#         Y = ETref_unstack.Y.values,
#         time_series_dates = ETref_unstack.time_series_dates.values,
#         realization = ETref_unstack.realization.values,
#     ),
#     attrs = dict(
#         Description = 'Reference crop evapotranspiration (mm/day)'),
# )

# var_OUT.ETo.values

# #Save as a netcdf for later processing
# var_OUT.to_netcdf(path = f'{home_dir}/ETo_lagged_average.nc4', mode ='w')
# #%%
# a = ETref_func.READ_DESCRIPTION_ATTRS[0,200,0].values
# print(a)


# #%%


# #%%

# for j in range(test.READ_DESCRIPTION_ATTRS.shape[1]):
#     for i in range(test.READ_DESCRIPTION_ATTRS.shape[2]):
#         if test.READ_DESCRIPTION_ATTRS[0,j,i] < 0:
#             print (test.READ_DESCRIPTION_ATTRS[0,i,i].values)
#             break
#         print (test.READ_DESCRIPTION_ATTRS[0,i,i].values) 
        
        
    



