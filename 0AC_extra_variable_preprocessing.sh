#!/bin/bash

main_directory='/home/kdl/Insync/OneDrive/NRT_CPC_Internship'
processors=8

data_d=$main_directory/Data
data_s=$main_directory/Scripts
mask=$data_s/CONUS_mask
subx=$data_d/SubX

model_array=(GMAO ESRL RSMAS EMC)

#Check latest dates downloaded
for mod in "${model_array[@]}";do
ls -ls $data_d/SubX/$mod | tail -7;done


############ Continue pre-processing SubX models ###############################
#Get the daily average, cdo operators do not do this properly because because the files
#were originally split. First 10 lead days. The 10th day had its data split between 2 seperate files.
python3 $data_s/EMC0_day_average_GEFSv12_HPC.py
python3 $data_s/EMC1_merge_3_soil_moisture_fields_GEFSv12.py

#Remove years 2011 and 2012 from the EMC GEFSv12 model (Deangelis et al. 2020)
cd $subx/EMC

rm *2011*
rm *2012*

#Priestley-Taylor coefficient -- COMPLETE
pt_dir=$data_d/Priestley_Taylor_Makkink_evap_coeff_maps
ncatted -O -a units,X,c,c,"degrees_east" -a units,Y,c,c,"degrees_south" pt_coeff.nc pt_coeff_add_metadata.nc

cdo -remapcon,$mask/CONUS_mask.grd $pt_dir/pt_coeff_add_metadata.nc $pt_dir/pt_coeff_FINAL.nc
################################################################################
#Elevation dataset preprocess (for altitude) -- COMPLETE
mkdir $data_d/elevation && cd $data_d/elevation
data_e=$data_d/elevation

cdo -remapcon,$mask/CONUS_mask.grd $data_e/elev.1-deg.nc $data_e/elev_regrid.nc
################################################################################

#Make evaporative demand datasets (All models)
python3 $data_s/ALL_MODEL_make_reference_ET.py



#Make evaporative demand dataset EMC
python3 $data_s/EMC3_make_anomaly_mean.py



#GMAO Soil moisture. Convert depth file given by Dr. Robert Koster to netcdf (python)
python3 $data_s/GMAO0_convert_depth_dat_file.py

#Remove unncessary dates GMAO
python3 $data_s/GMAO1_remove_unncessary_dates.py

#Convert soil moisture GMAO to proper format (m3/m3)
python3 $data_s/GMAO2_convert_RZSM_to_m3_m3.py

#Make evaporative demand dataset GMAO
python3 $data_s/GMAO3_make_reference_ET.py

#Create a mean dataset for anomaly
python3 $data_s/GMAO4_make_anomaly_mean.py



#Make empty anomaly mean files to store output




